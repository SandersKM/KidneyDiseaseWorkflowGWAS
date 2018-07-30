
start_time <- Sys.time()
# source('http://bioconductor.org/biocLite.R') biocLite('GenomeGraphs') install.packages('devtools')
# devtools::install_github('rlbarter/superheat') install_github('arcdiagram', username='gastonstat')
library(jsonlite)
library(xml2)
library(rentrez)
library(gwascat)
library(magrittr)
library(GenomeGraphs)
library(biomaRt)
library(ensembldb)
library(ggplot2)
library(httr)
library(ComplexHeatmap)
library(GetoptLong)
library(circlize)
library(stringi)
library(data.table)
library(arcdiagram)
library(scales)

###################
# User Input
###################

# Enter the HGNC symbol of the gene of interest below
gene.of.interest <- "MUC1"

# Enter the filepath you would like all documents produced to go to
output.dir <- "/Users/ksanders/Documents/"

# This will bring you to the nephvs eQTL webpage
browseURL(paste("http://eqtl.nephvs.org/searchResult/", gene.of.interest, sep = ""))
# Enter the path to the location of these downloaded files
filepath <- "/Users/ksanders/Downloads/"

# Read in the NephQTL data Change filename (arg 2) if it is different from below
nephQTL.glom <- read.csv(paste(filepath, "glom_MatrixEQTL_", gene.of.interest, ".csv", sep = ""), header = TRUE, sep = ",")
nephQTL.tub <- read.csv(paste(filepath, "tub_MatrixEQTL_", gene.of.interest, ".csv", sep = ""), header = TRUE, sep = ",")

# 1000 Genomes Project population for LD retrieval
# Descriptions of populations: http://grch37.rest.ensembl.org/documentation/info/variation_populations
# Currently using Utah 1000 genomes data
genomes.population <- "CEU"

# Chose ONE way to sort the output file
sort.by.pvalue <- TRUE # sorted by highest significance accross all tissue types
sort.by.GWAS <- FALSE # variants with the most GWAS hits are first.
sort.by.position <- FALSE # sorted by bp order

####################
# BioMart Gene Info
####################

# Setting up biomart
if (!exists("mart")) {
    mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37)
}
gene.of.interest.info <- getBM(c("start_position", "end_position", "strand", "chromosome_name"), filters = "hgnc_symbol", values = gene.of.interest,
    mart = mart)
gene.of.interest.start <- gene.of.interest.info$start_position
gene.of.interest.strand <- gene.of.interest.info$strand
gene.of.interest.end <- gene.of.interest.info$end_position
gene.of.interest.chrom <- gene.of.interest.info$chromosome_name

########################################
# Make table of eQTL positions/values
########################################

# Initialize combined dataframe
total.rows <- dim(nephQTL.glom)[1] + dim(nephQTL.tub)[1]
eQTL.combined <- data.frame(SNPid = character(total.rows), chrom = character(total.rows), position = character(total.rows), ref = character(total.rows),
    alt = character(total.rows), pvalue = character(total.rows), beta = character(total.rows), compartment = character(total.rows), source = character(total.rows),
    stringsAsFactors = FALSE)

# filling with nephQTL glom and tub tables
sapply(1:dim(nephQTL.glom)[1], function(n) {
    eQTL.combined$SNPid[n] <<- toString(nephQTL.glom$dbSNPId[n])
    chr_pos <- strsplit(toString(nephQTL.glom$Chr.pos[n]), split = ":")
    eQTL.combined$chrom[n] <<- chr_pos[[1]][1]
    eQTL.combined$position[n] <<- chr_pos[[1]][2]
    eQTL.combined$ref[n] <<- toString(nephQTL.glom$Ref.[n])
    eQTL.combined$alt[n] <<- toString(nephQTL.glom$Alt.[n])
    eQTL.combined$pvalue[n] <<- toString(nephQTL.glom$P.value[n])
    eQTL.combined$beta[n] <<- toString(nephQTL.glom$Beta[n])
    eQTL.combined$compartment[n] <<- "Glom"
    eQTL.combined$source[n] <<- "NephQTL"
})
sapply(1:dim(nephQTL.tub)[1], rowstart = dim(nephQTL.glom)[1], function(n, rowstart) {
    i <- rowstart + n
    eQTL.combined$SNPid[i] <<- toString(nephQTL.tub$dbSNPId[n])
    chr_pos <- strsplit(toString(nephQTL.tub$Chr.pos[n]), split = ":")
    eQTL.combined$chrom[i] <<- chr_pos[[1]][1]
    eQTL.combined$position[i] <<- chr_pos[[1]][2]
    eQTL.combined$ref[i] <<- toString(nephQTL.tub$Ref.[n])
    eQTL.combined$alt[i] <<- toString(nephQTL.tub$Alt.[n])
    eQTL.combined$pvalue[i] <<- toString(nephQTL.tub$P.value[n])
    eQTL.combined$beta[i] <<- toString(nephQTL.tub$Beta[n])
    eQTL.combined$compartment[i] <<- "Tub"
    eQTL.combined$source[i] <<- "NephQTL"
})

eqtl.combined.tub <- eQTL.combined[which(eQTL.combined$compartment == "Tub"),]
eqtl.combined.tub$pvalue <- as.numeric(as.character(eqtl.combined.tub$pvalue))
eqtl.combined.tub$Mlog <- -log10(eqtl.combined.tub$pvalue)
eqtl.combined.tub.length <- dim(eqtl.combined.tub)[1]
eqtl.combined.glom <- eQTL.combined[which(eQTL.combined$compartment == "Glom"),]
eqtl.combined.glom$pvalue <- as.numeric(as.character(eqtl.combined.glom$pvalue))
eqtl.combined.glom$Mlog <- -log10(eqtl.combined.glom$pvalue)
eqtl.combined.glom.length <- dim(eqtl.combined.glom)[1]


# make Mlog for eQTL.combined

eQTL.combined$pvalue <- as.numeric(as.character(eQTL.combined$pvalue))
eQTL.combined$Mlog <- -log10(eQTL.combined$pvalue)

######################################################################
# Get GWAS results for variants reported or mapped to Gene of Interest
######################################################################

# Loading current GWAS Catelogue information.  Do this every few months, but be warned it takes a long time to load
if (!exists("current_gwascat")) {
    current_gwascat <- makeCurrentGwascat(genome = "GRCh37")
}

get_reported_gene_of_interest <- function(numreported) {
    rows <- ""
    for (i in 1:numreported) {
        isReported <- gene.of.interest %in% strsplit(gsub(" ", "", current_gwascat$`REPORTED GENE(S)`[i]), split = ",")[[1]]
        if (isReported) {
            rows <- paste(rows, i, ",", sep = "")
        }
    }
    return(strsplit(rows, split = ","))
}

reported_gene_rows <- as.integer(get_reported_gene_of_interest(length(current_gwascat$`REPORTED GENE(S)`))[[1]])
mapped_gene_rows <- which(current_gwascat$MAPPED_GENE == gene.of.interest)
combined_gene_rows <- unique(append(mapped_gene_rows, reported_gene_rows))

# Data frame with all gwas hits mapped or reported to gene of interest
gwas.all.hits <- as.data.frame(current_gwascat[combined_gene_rows])

# Condensed version of gwas.all.hits, where each variant postion has a row
gwas.variants <- data.frame(RSID = unique(gwas.all.hits$SNPS))
gwas.variants$chrom <- character(dim(gwas.variants)[1])
gwas.variants$position <- character(dim(gwas.variants)[1])
gwas.variants$all.studies <- character(dim(gwas.variants)[1])
gwas.variants$risk.allele.A <- character(dim(gwas.variants)[1])
gwas.variants$risk.allele.C <- character(dim(gwas.variants)[1])
gwas.variants$risk.allele.G <- character(dim(gwas.variants)[1])
gwas.variants$risk.allele.T <- character(dim(gwas.variants)[1])
gwas.variants$risk.allele.NA <- character(dim(gwas.variants)[1])
gwas.variants$pubmed.links <- character(dim(gwas.variants)[1])
get_risk_allele_string <- function(r, a) {
    allele.rows <- gwas.all.hits.rows[which(substring(r$STRONGEST.SNP.RISK.ALLELE, first = stri_length(r$STRONGEST.SNP.RISK.ALLELE[1])) ==
        a), ]
    return(paste("Disease Trait: ", allele.rows$DISEASE.TRAIT, "; P-Value: ", allele.rows$P.VALUE, "; -log_10(P-Value): ", allele.rows$PVALUE_MLOG,
        "; Odds Ratio / BETA: ", allele.rows$OR.or.BETA, "; 95% CI: ", allele.rows$X95..CI..TEXT., collapse = " || "))
}
for (i in 1:dim(gwas.variants)[1]) {
    gwas.all.hits.rows <- gwas.all.hits[which(gwas.all.hits$SNPS == gwas.variants$RSID[i]), ]
    gwas.variants$chrom[i] <- gwas.all.hits.rows$CHR_ID[1]
    gwas.variants$position[i] <- gwas.all.hits.rows$start[1]
    gwas.variants$all.studies[i] <- paste(unique(gwas.all.hits.rows$STUDY), collapse = ";    ")
    gwas.variants$risk.allele.A[i] <- get_risk_allele_string(gwas.all.hits.rows, "A")
    gwas.variants$risk.allele.C[i] <- get_risk_allele_string(gwas.all.hits.rows, "C")
    gwas.variants$risk.allele.G[i] <- get_risk_allele_string(gwas.all.hits.rows, "G")
    gwas.variants$risk.allele.T[i] <- get_risk_allele_string(gwas.all.hits.rows, "T")
    gwas.variants$risk.allele.NA[i] <- get_risk_allele_string(gwas.all.hits.rows, "?")
    gwas.variants$pubmed.links[i] <- paste(gwas.all.hits.rows$LINK, collapse = " ; ")
}

##################
# Getting LD Data
##################

# Variants that are GWAS hits with significant eQTL
eQTL.gwas.combined.rsid <- intersect(gwas.variants$RSID, eQTL.combined$SNPid)
eQTL.gwas.combined.LD.r2 <- matrix(nrow = length(eQTL.gwas.combined.rsid), ncol = length(eQTL.gwas.combined.rsid))
eQTL.gwas.combined.LD.dprime <- matrix(nrow = length(eQTL.gwas.combined.rsid), ncol = length(eQTL.gwas.combined.rsid))
colnames(eQTL.gwas.combined.LD.r2) <- colnames(eQTL.gwas.combined.LD.dprime) <- eQTL.gwas.combined.rsid
rownames(eQTL.gwas.combined.LD.r2) <- rownames(eQTL.gwas.combined.LD.dprime) <- eQTL.gwas.combined.rsid
for (i in 1:(length(eQTL.gwas.combined.rsid) - 1)) {
    for (j in (i + 1):length(eQTL.gwas.combined.rsid)) {
        ensembl.ld.json <- read_json(paste("http://grch37.rest.ensembl.org/ld/human/pairwise/", eQTL.gwas.combined.rsid[i], "/", eQTL.gwas.combined.rsid[j],
            "?content-type=application/json;population_name=1000GENOMES:", "phase_3:", genomes.population, sep = ""))
        if (length(ensembl.ld.json) > 0) {
            eQTL.gwas.combined.LD.r2[i, j] <- eQTL.gwas.combined.LD.r2[j, i] <- ensembl.ld.json[[1]]$r2
            eQTL.gwas.combined.LD.dprime[i, j] <- eQTL.gwas.combined.LD.dprime[j, i] <- ensembl.ld.json[[1]]$d_prime
        }
    }
}

# LD data within 500 KB of the center of the gene of interest
gene.of.interest.half <- (gene.of.interest.start + gene.of.interest.end)%/%2
start.500Kb <- gene.of.interest.half - 250000
end.500Kb <- gene.of.interest.half + 250000
LD.info.500Kb <- read_json(paste("http://grch37.rest.ensembl.org/ld/human/region/", gene.of.interest.chrom, ":", start.500Kb, "..", end.500Kb -
    1, "/1000GENOMES:phase_3:", genomes.population, "?content-type=application/json", sep = ""))
# make nested list into data frame
LD.info.500Kb <- rbindlist(LD.info.500Kb, fill = FALSE)

num.unique <- length(union(LD.info.500Kb$variation1, LD.info.500Kb$variation2))

LD.info.500Kb.unique.variants <- data.frame(rsid = union(LD.info.500Kb$variation1, LD.info.500Kb$variation2),
                                            chrom = gene.of.interest.chrom,
                                            position = numeric(num.unique), stringsAsFactors = FALSE)


# Fetch Position from Ensembl using RSID
server <- "http://grch37.rest.ensembl.org"
ext <- "/vep/human/id"
i <- 1
while (i < dim(LD.info.500Kb.unique.variants)[1]) {
    j <- i + 190  # Ensembl takes at most 200 requests at a time.
    if (j > dim(LD.info.500Kb.unique.variants)[1]) {
        j = dim(LD.info.500Kb.unique.variants)[1]
    }
    rest.api.response <- POST(paste(server, ext, sep = ""), content_type("application/json"), accept("application/json"), body = paste("{ \"ids\" : [",
        paste0(LD.info.500Kb.unique.variants$rsid[i:j], collapse = "\",\""), " ] }", sep = "\""))
    rest.api.info <- as.data.frame(fromJSON(toJSON(content(rest.api.response))))
    for (k in 1:dim(rest.api.info)[1]) {
        if (length(rest.api.info$end[k][[1]]) > 0) {
            LD.info.500Kb.unique.variants$position[k + i - 1] <- unlist(rest.api.info$end[k])
        } else {
            LD.info.500Kb.unique.variants$position[k + i - 1] <- 0
        }
    }
    i <- j + 1
}

# Get rid of 0s if applicable and sort ascending by position
LD.info.500Kb.unique.variants <- LD.info.500Kb.unique.variants[which(LD.info.500Kb.unique.variants$position > 0), ]
LD.info.500Kb.unique.variants <- LD.info.500Kb.unique.variants[order(LD.info.500Kb.unique.variants$position), ]

# full dataframe where LD and eQTL overlap
ld.eqtl.overlap <- LD.info.500Kb[which(LD.info.500Kb$variation1 %in% eQTL.combined$SNPid &
                                         LD.info.500Kb$variation2 %in% eQTL.combined$SNPid),]
ld.eqtl.overlap$position.variation1 <- LD.info.500Kb.unique.variants$position[match(ld.eqtl.overlap$variation1, LD.info.500Kb.unique.variants$rsid)]
ld.eqtl.overlap$position.variation2 <- LD.info.500Kb.unique.variants$position[match(ld.eqtl.overlap$variation2, LD.info.500Kb.unique.variants$rsid)]

ld.eqtl.overlap$eQTL.Tub.pvalue.variation2 <- eqtl.combined.tub$pvalue[match(ld.eqtl.overlap$variation2, eqtl.combined.tub$SNPid)]
ld.eqtl.overlap$eQTL.Tub.pvalue.variation1 <- eqtl.combined.tub$pvalue[match(ld.eqtl.overlap$variation1, eqtl.combined.tub$SNPid)]
ld.eqtl.overlap$eQTL.Glom.pvalue.variation2 <- eqtl.combined.glom$pvalue[match(ld.eqtl.overlap$variation2, eqtl.combined.glom$SNPid)]
ld.eqtl.overlap$eQTL.Glom.pvalue.variation1 <- eqtl.combined.glom$pvalue[match(ld.eqtl.overlap$variation1, eqtl.combined.glom$SNPid)]

####################
# Make output file
####################

# Starting with all unique positions in eQTL data
unique.positions <- unique(eQTL.combined$position)
unique.position.rownum <- match(unique.positions,eQTL.combined$position)
output.file <- data.frame(rsid = eQTL.combined$SNPid[unique.position.rownum],
                          chrom = gene.of.interest.chrom, position = unique.positions,
                          ref = eQTL.combined$ref[unique.position.rownum],
                          alt = eQTL.combined$alt[unique.position.rownum],
                          stringsAsFactors = FALSE)
eqtl.combined.glom.row <- match(output.file$rsid,eqtl.combined.glom$SNPid)
eqtl.combined.tub.row <- match(output.file$rsid,eqtl.combined.tub$SNPid)
output.file$tub.pvalue <- eqtl.combined.tub$pvalue[eqtl.combined.tub.row]
output.file$tub.MLog <- eqtl.combined.tub$Mlog[eqtl.combined.tub.row]
output.file$tub.beta <- eqtl.combined.tub$beta[eqtl.combined.tub.row]
output.file$glom.pvalue <- eqtl.combined.glom$pvalue[eqtl.combined.glom.row]
output.file$glom.MLog <- eqtl.combined.glom$Mlog[eqtl.combined.glom.row]
output.file$glom.beta <- eqtl.combined.glom$beta[eqtl.combined.glom.row]

# To get a row that has GWAS information for that variant

output.file$GWAS <- character(length(output.file$rsid))
gwas.match.rownum <- match(gwas.variants$RSID, output.file$rsid)
for(i in 1:length(gwas.variants$RSID)){
  if(!is.na(gwas.match.rownum[i])){
    output.file$GWAS[gwas.match.rownum[i]] <- gwas.variants$all.studies[i]
  }
}

# To get a row that has all LD informationrelated to that variant
add_ld_information_to_output <- function(rsid){
  ld.string <- ""
  var1.rownum <- which(rsid==LD.info.500Kb$variation1)
  var2.rownum <- which(rsid==LD.info.500Kb$variation2)
  if(length(union(var1.rownum, var2.rownum)) ==0){
    return("")
  }
  ld.string <- paste0(LD.info.500Kb$variation2[which(rsid==LD.info.500Kb$variation1)], " : R2 = ",
                      LD.info.500Kb$r2[ which(rsid==LD.info.500Kb$variation1)], " : D_prime = ",
                      LD.info.500Kb$d_prime[ which(rsid==LD.info.500Kb$variation1)],"; ", collapse = "")
  ld.string <- paste(ld.string, paste0(LD.info.500Kb$variation1[which(rsid==LD.info.500Kb$variation2)], " : R2 = ",
                                        LD.info.500Kb$r2[ which(rsid==LD.info.500Kb$variation2)], " : D_prime = ",
                                        LD.info.500Kb$d_prime[ which(rsid==LD.info.500Kb$variation2)],"; ", collapse = ""),sep = "")
  return(ld.string)
}
output.file$LD.info <- sapply(output.file$rsid, add_ld_information_to_output)

# Sorting output file

if(sort.by.position){
  output.file <- output.file[order(output.file$position),]
}

if(sort.by.GWAS){
  output.file$temp <- sapply(output.file$GWAS,function(x){length(strsplit(x, split = ";")[[1]])})
  output.file <- output.file[order(output.file$temp, decreasing = TRUE),]
  output.file <- output.file[, !colnames(output.file) == "temp"]
}

if(sort.by.pvalue){
  output.file$temp <- sapply(1:dim(output.file)[1], function(n){
    return(max(output.file$tub.MLog[n], output.file$glom.MLog[n]))
  })
  output.file <- output.file[order(output.file$temp, decreasing = TRUE),]
  output.file <- output.file[, !colnames(output.file) == "temp"]
}

#####################
# Gene Plot all eQTL
#####################

minbase <- min(as.numeric(eQTL.combined$position))
maxbase <- max(as.numeric(eQTL.combined$position))
gene.image <- makeGene(id = gene.of.interest, type = "hgnc_symbol", biomart = mart)
genesplus <- makeGeneRegion(start = minbase, end = maxbase, strand = "+", chromosome = gene.of.interest.chrom, biomart = mart)
genesmin <- makeGeneRegion(start = minbase, end = maxbase, strand = "-", chromosome = gene.of.interest.chrom, biomart = mart)
expres.glom <- makeSegmentation(value = as.numeric(-log10(as.numeric(eQTL.combined$pvalue[which(eQTL.combined$compartment == "Glom")]))),
    start = as.numeric(eQTL.combined$position[which(eQTL.combined$compartment == "Glom")]), end = as.numeric(eQTL.combined$position[which(eQTL.combined$compartment ==
        "Glom")]), dp = DisplayPars(color = "#f2b229", lwd = 7, lty = 1))
expres.tub <- makeGenericArray(intensity = as.matrix(-log10(as.numeric(eQTL.combined$pvalue[which(eQTL.combined$compartment == "Tub")]))),
    probeStart = as.numeric(eQTL.combined$position[which(eQTL.combined$compartment == "Tub")]), dp = DisplayPars(type = "point", pch = 19,
        color = "#a162c4"), trackOverlay = expres.glom)
setPar(expres.tub, "pointSize", 0.7)
genomeAxis <- makeGenomeAxis(add53 = TRUE, add35 = TRUE)
gene.region.overlay <- makeRectangleOverlay(start = gene.of.interest.start, end = gene.of.interest.end, dp = DisplayPars(fill = "yellow",
    alpha = 0.2, lty = "dashed"))
legend = makeLegend(text = c("Tub", "Glom", "GWAS", gene.of.interest), fill = c("#a162c4", "#f2b229", "#2872f1", "lightyellow"), cex = 1)

overlays <- vector(mode = "list", length = dim(gwas.variants)[1] + 1)
overlays[dim(gwas.variants)[1] + 1] <- makeRectangleOverlay(start = gene.of.interest.start, end = gene.of.interest.end, dp = DisplayPars(fill = "yellow",
    alpha = 0.2, lty = "dotted"), region = c(2, 3))
for (i in 1:dim(gwas.variants)[1]) {
    overlays[i] <- makeTextOverlay("o", xpos = as.numeric(gwas.variants$position[i]), ypos = 0.075, coords = "genomic", dp = DisplayPars(color = "#2872f1",
        cex = 2.5))
}
jpeg(file = paste(output.dir, gene.of.interest, "-eQTL_GWAS", ".jpeg", sep = ""), width = 1000, height = 900)
gdPlot(list(legend, `-log(P value)  of  eQTL` = expres.tub, BP = genomeAxis), overlays = overlays, minBase = minbase, maxBase = maxbase,
    labelCex = 2)
dev.off()

#######################################
# LD Heatmap and eQTL around each GWAS
#######################################


ld.eqtl.gwas.overlap.position <- LD.info.500Kb.unique.variants[match(eQTL.gwas.combined.rsid, LD.info.500Kb.unique.variants$rsid), ]

range.around.gwas <- 10000
dprime.min <- 0.85


for (n in 1:length(ld.eqtl.gwas.overlap.position$position)) {
  if(!is.na(ld.eqtl.gwas.overlap.position$position[n])){
    zoom.ld.eqtl.gwas.overlap <- ld.eqtl.overlap[which(ld.eqtl.overlap$position.variation1 < ld.eqtl.gwas.overlap.position$position[n] +
        range.around.gwas & ld.eqtl.overlap$position.variation1 > ld.eqtl.gwas.overlap.position$position[n] - range.around.gwas & ld.eqtl.overlap$position.variation2 <
        ld.eqtl.gwas.overlap.position$position[n] + range.around.gwas & ld.eqtl.overlap$position.variation2 > ld.eqtl.gwas.overlap.position$position[n] -
        range.around.gwas), ]

    zoom.ld.eqtl.gwas.overlap <- zoom.ld.eqtl.gwas.overlap[(zoom.ld.eqtl.gwas.overlap$d_prime > dprime.min | zoom.ld.eqtl.gwas.overlap$variation1 %in%
        ld.eqtl.gwas.overlap.position$rsid | zoom.ld.eqtl.gwas.overlap$variation2 %in% ld.eqtl.gwas.overlap.position$rsid), ]

    # number of LD - eQTL overlaps in zoom
    zoom.ld.eqtl.gwas.overlap.row.num <- length(unique(union(zoom.ld.eqtl.gwas.overlap$variation1, zoom.ld.eqtl.gwas.overlap$variation2)))
    # rsids of overlap in zoom
    zoom.ld.eqtl.gwas.overlap.rsids <- unique(union(zoom.ld.eqtl.gwas.overlap$variation1, zoom.ld.eqtl.gwas.overlap$variation2))

    # sorted rsids of overlap
    zoom.ld.eqtl.gwas.overlap.rsids <- zoom.ld.eqtl.gwas.overlap.rsids[order(LD.info.500Kb.unique.variants$position[match(zoom.ld.eqtl.gwas.overlap.rsids,
        LD.info.500Kb.unique.variants$rsid)])]

    # BP postions in LD - eQTL overlap sorted smallest to largest
    zoom.ld.eqtl.gwas.overlap.positions <- LD.info.500Kb.unique.variants$position[match(zoom.ld.eqtl.gwas.overlap.rsids, LD.info.500Kb.unique.variants$rsid)]
    zoom.ld.eqtl.gwas.overlap.total.dist <- max(zoom.ld.eqtl.gwas.overlap.positions) - min(zoom.ld.eqtl.gwas.overlap.positions)

    # matrix for distances
    zoom.ld.eqtl.gwas.overlap.dist <- outer(zoom.ld.eqtl.gwas.overlap.positions, zoom.ld.eqtl.gwas.overlap.positions, "-")
    zoom.ld.eqtl.gwas.overlap.dist[upper.tri(zoom.ld.eqtl.gwas.overlap.dist)] = 0
    zoom.ld.eqtl.gwas.overlap.dist.scaled <- rescale(zoom.ld.eqtl.gwas.overlap.dist, to = c(0, 1))

    # initialize matrix for r2 values
    zoom.ld.eqtl.gwas.overlap.r2 <- matrix(0, nrow = zoom.ld.eqtl.gwas.overlap.row.num, ncol = zoom.ld.eqtl.gwas.overlap.row.num, dimnames = list(zoom.ld.eqtl.gwas.overlap.rsids,
        zoom.ld.eqtl.gwas.overlap.rsids))

    # populate the r2 matrix
    for (i in 1:dim(zoom.ld.eqtl.gwas.overlap)[1]) {
        zoom.ld.eqtl.gwas.overlap.r2[zoom.ld.eqtl.gwas.overlap$variation1[i], zoom.ld.eqtl.gwas.overlap$variation2[i]] <- zoom.ld.eqtl.gwas.overlap.r2[zoom.ld.eqtl.gwas.overlap$variation2[i],
            zoom.ld.eqtl.gwas.overlap$variation1[i]] <- zoom.ld.eqtl.gwas.overlap$r2[i]
    }

    # making lower portion of LD plot white
    zoom.ld.eqtl.gwas.overlap.r2[lower.tri(zoom.ld.eqtl.gwas.overlap.r2)] = 0

    # For the top of the heatmap. This dataframe has the Tub and Glom eQTL data
    annot.df <- as.data.frame(cbind(eqtl.combined.tub$Mlog[match(rownames(zoom.ld.eqtl.gwas.overlap.r2), eqtl.combined.tub$SNPid)], eqtl.combined.glom$Mlog[match(rownames(zoom.ld.eqtl.gwas.overlap.r2),
        eqtl.combined.glom$SNPid)]), stringsAsFactors = FALSE)
    names(annot.df) <- c("Tubulointerstitium", "Glomerulus")
    annot.df$Tubulointerstitium <- as.numeric(annot.df$Tubulointerstitium)
    annot.df$Glomerulus <- as.numeric(annot.df$Glomerulus)
    annot <- HeatmapAnnotation(df = annot.df, col = list(Tubulointerstitium = colorRamp2(c(0, 6), c("white", "darkblue")), Glomerulus = colorRamp2(c(0,
        6), c("white", "orange"))), show_annotation_name = TRUE, show_legend = FALSE)


    # Get matrix with the combined values for the plot
    zoom.ld.eqtl.gwas.overlap.combined.r2.dist <- zoom.ld.eqtl.gwas.overlap.r2
    zoom.ld.eqtl.gwas.overlap.combined.r2.dist[lower.tri(zoom.ld.eqtl.gwas.overlap.combined.r2.dist)] <- zoom.ld.eqtl.gwas.overlap.dist.scaled[lower.tri(zoom.ld.eqtl.gwas.overlap.dist.scaled)]
    class(zoom.ld.eqtl.gwas.overlap.combined.r2.dist) <- "numeric"

    # color scale for r2 and distances
    col_r2 = colorRamp2(c(0, 1), c("white", "darkred"))
    col_dist = colorRamp2(c(0, 1), c("white", "black"))

    col_label_names <- match(rownames(zoom.ld.eqtl.gwas.overlap.combined.r2.dist), as.character(ld.eqtl.gwas.overlap.position$rsid))
    col_label_names[!is.na(col_label_names)] <- "#2872f1"
    col_label_names[is.na(col_label_names)] <- "#404040"

    # legend information
    legend_Tubulointerstitium = Legend(at = seq(0, max(annot.df), by = 1), title = "Tubulointerstitium\n-log(pvalue)", col_fun = colorRamp2(c(0,
        6), c("white", "darkblue")), legend_height = unit(1, "cm"))
    legend_Glomerulus = Legend(at = seq(0, max(annot.df), by = 1), col_fun = colorRamp2(c(0, 6), c("white", "orange")), title = "\nGlomerulus\n-log(pvalue)")
    legend_r2 = Legend(at = seq(0, 1, by = 0.2), title = "\nLinkage\nDisequilibrium\nR2", col_fun = colorRamp2(c(0, 1), c("white", "darkred")))
    legend_dist = Legend(at = seq(0, max(zoom.ld.eqtl.gwas.overlap.dist), by = 10000), title = "\nDistance", col_fun = colorRamp2(c(0,  max(zoom.ld.eqtl.gwas.overlap.dist)), c("white", "black")))
    legend_gwas = Legend(at = c("GWAS"), title = "GWAS HIT", type = "points", legend_gp = gpar(col = "#2872f1"))

    # Make the heatmap!
    zoom.ld.eqtl.gwas.overlap.Heatmap <- Heatmap(zoom.ld.eqtl.gwas.overlap.combined.r2.dist, col = colorRamp2(c(-1, 1), c("white", "white"),
        transparency = 0.5), name = "LD R2", cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE, show_row_names = TRUE,
        row_names_gp = gpar(col = col_label_names), show_column_names = TRUE, column_names_gp = gpar(col = col_label_names), show_column_dend = TRUE,
        top_annotation = annot, top_annotation_height = unit(4, "cm"), cell_fun = function(j, i, x, y, width, height, fill) {
            grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "white", fill = NA))
            if (i < j) {
                grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "lightgrey", fill = col_r2(zoom.ld.eqtl.gwas.overlap.combined.r2.dist[i,
                  j])))
            } else if (i > j) {
                grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = col_dist(zoom.ld.eqtl.gwas.overlap.combined.r2.dist[i,
                  j]), fill = col_dist(zoom.ld.eqtl.gwas.overlap.combined.r2.dist[i, j])))
            }
        })
    jpeg(file = paste(output.dir, gene.of.interest, "-GWAS-", ld.eqtl.gwas.overlap.position$rsid[n], "-LD_eQTL", ".jpeg", sep = ""), width = 1000,
        height = 900)
    draw(zoom.ld.eqtl.gwas.overlap.Heatmap, annotation_legend_list = list(legend_Tubulointerstitium, legend_Glomerulus, legend_r2, legend_dist,
        legend_gwas))
    dev.off()
  }
}
end_time <- Sys.time()

####################
# Write Output File
####################

write.csv(output.file, file = paste(output.dir, gene.of.interest, "_ComplexVariants.csv", sep = ""), row.names = FALSE)



