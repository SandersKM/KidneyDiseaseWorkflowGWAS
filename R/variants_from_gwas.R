###################
# Kate Sanders
# BSRP 2018
###################

# source("http://bioconductor.org/biocLite.R")
# biocLite("GenomeGraphs")
# install.packages("devtools")
# devtools::install_github("rlbarter/superheat")
# install_github('arcdiagram',  username='gastonstat')
library(jsonlite)
library(rentrez)
library(gwascat)
library(magrittr)
library(GenomeGraphs)
library(ensembldb)
library(ggplot2)
library(httr)
library(ComplexHeatmap)
library(circlize)
library(stringi)
library(arcdiagram)

# Enter the HGNC symbol of the gene of interest below
gene.of.interest <- "MUC1"

# Enter the filepath you would like all documents produced to go to
gwas.filepath <- "/Users/ksanders/Documents/"

# This will bring you to the nephvs eQTL webpage
browseURL(paste("http://eqtl.nephvs.org/searchResult/", gene.of.interest, sep = ""))
# Enter the path to the location of these downloaded files
filepath <- "/Users/ksanders/Desktop/"

# Read in the NephQTL data
nephQTL.glom <- read.csv(paste(filepath, "glom_MatrixEQTL_", gene.of.interest,".csv", sep = ""),
                         header = TRUE, sep = ",")
nephQTL.tub <- read.csv(paste(filepath, "tub_MatrixEQTL_", gene.of.interest,".csv", sep = ""),
                        header = TRUE, sep = ",")
# Setting up biomart
if(!exists("mart")){
  mart <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
}
gene.of.interest.info <- getBM(c("start_position", "end_position", "strand", "chromosome_name"),
                               filters="hgnc_symbol",values=gene.of.interest, mart=mart)
gene.of.interest.start <- gene.of.interest.info$start_position
gene.of.interest.strand <- gene.of.interest.info$strand
gene.of.interest.end <- gene.of.interest.info$end_position
gene.of.interest.chrom <- gene.of.interest.info$chromosome_name

##############################################
# Make table of eQTL positions/values
##############################################

total.rows <- dim(nephQTL.glom)[1] + dim(nephQTL.tub)[1]
eQTL.combined <- data.frame(SNPid = character(total.rows), chrom = character(total.rows),
                            position = character(total.rows), ref = character(total.rows),
                            alt = character(total.rows), pvalue = character(total.rows),
                            beta = character(total.rows), compartment = character(total.rows),
                            source = character(total.rows), stringsAsFactors = FALSE)

# filling with nephQTL glom and tub tables
sapply(1:dim(nephQTL.glom)[1], function(n){
  eQTL.combined$SNPid[n] <<- toString(nephQTL.glom$dbSNPId[n])
  chr_pos <- strsplit(toString(nephQTL.glom$Chr.pos[n]), split = ":")
  eQTL.combined$chrom[n]<<-chr_pos[[1]][1]
  eQTL.combined$position[n]<<-chr_pos[[1]][2]
  eQTL.combined$ref[n]<<- toString(nephQTL.glom$Ref.[n])
  eQTL.combined$alt[n]<<-toString(nephQTL.glom$Alt.[n])
  eQTL.combined$pvalue[n]<<-toString(nephQTL.glom$P.value[n])
  eQTL.combined$beta[n]<<-toString(nephQTL.glom$Beta[n])
  eQTL.combined$compartment[n]<<-"Glom"
  eQTL.combined$source[n] <<-"NephQTL"
})
sapply(1:dim(nephQTL.tub)[1], rowstart = dim(nephQTL.glom)[1], function(n, rowstart){
  i <- rowstart + n
  eQTL.combined$SNPid[i] <<- toString(nephQTL.tub$dbSNPId[n])
  chr_pos <- strsplit(toString(nephQTL.tub$Chr.pos[n]), split = ":")
  eQTL.combined$chrom[i]<<-chr_pos[[1]][1]
  eQTL.combined$position[i]<<-chr_pos[[1]][2]
  eQTL.combined$ref[i]<<- toString(nephQTL.tub$Ref.[n])
  eQTL.combined$alt[i]<<-toString(nephQTL.tub$Alt.[n])
  eQTL.combined$pvalue[i]<<-toString(nephQTL.tub$P.value[n])
  eQTL.combined$beta[i]<<-toString(nephQTL.tub$Beta[n])
  eQTL.combined$compartment[i]<<-"Tub"
  eQTL.combined$source[i] <<-"NephQTL"
})

######################################################################
# Get GWAS results for variants reported or mapped to Gene of Interest
######################################################################

# Loading current GWAS Cateloge information.
# Do this every few months, but be warned it take a long time to load
if(!exists("current_gwascat")){
  current_gwascat <- makeCurrentGwascat(genome = "GRCh37")
}

get_reported_gene_of_interest <- function(numreported){
  rows <-  ""
  for(i in 1:numreported){
    isReported <- gene.of.interest %in%  strsplit(gsub(" ", "", current_gwascat$`REPORTED GENE(S)`[i])
                                                  , split = ",")[[1]]
    if(isReported){
      rows <- paste(rows , i , ",",sep = "")
    }
  }
  return(strsplit(rows, split = ","))
}

reported_gene_rows <- as.integer(get_reported_gene_of_interest(length(
  current_gwascat$`REPORTED GENE(S)`))[[1]])
mapped_gene_rows <- which(current_gwascat$MAPPED_GENE == gene.of.interest)
combined_gene_rows <- unique(append(mapped_gene_rows, reported_gene_rows))
# Data frame with all gwas hits mapped or reported to gene of interest
gwas.all.hits <- as.data.frame(current_gwascat[combined_gene_rows])

# Condensed version of gwas.all.hits, where each variant postion has a row
gwas.variants <- data.frame(RSID = unique(gwas.all.hits$SNPS))
gwas.variants$chrom <- character(dim(gwas.variants)[1])
gwas.variants$position <- character(dim(gwas.variants)[1])
gwas.variants$all.studies  <- character(dim(gwas.variants)[1])
gwas.variants$risk.allele.A <- character(dim(gwas.variants)[1])
gwas.variants$risk.allele.C  <- character(dim(gwas.variants)[1])
gwas.variants$risk.allele.G  <- character(dim(gwas.variants)[1])
gwas.variants$risk.allele.T <- character(dim(gwas.variants)[1])
gwas.variants$risk.allele.NA  <- character(dim(gwas.variants)[1])
gwas.variants$pubmed.links  <- character(dim(gwas.variants)[1])
get_risk_allele_string <- function(r, a){
  allele.rows <- gwas.all.hits.rows[which(substring(
    r$STRONGEST.SNP.RISK.ALLELE, first = stri_length(r$STRONGEST.SNP.RISK.ALLELE[1])) == a),]
  return(paste("Disease Trait: ", allele.rows$DISEASE.TRAIT,"; P-Value: ", allele.rows$P.VALUE,
               "; -log_10(P-Value): ", allele.rows$PVALUE_MLOG, "; Odds Ratio / BETA: ", allele.rows$OR.or.BETA,
               "; 95% CI: ", allele.rows$X95..CI..TEXT., collapse = " | "))
}
for(i in 1:dim(gwas.variants)[1]){
  gwas.all.hits.rows <- gwas.all.hits[which(gwas.all.hits$SNPS == gwas.variants$RSID[i]),]
  gwas.variants$chrom[i] <- gwas.all.hits.rows$CHR_ID[1]
  gwas.variants$position[i] <- gwas.all.hits.rows$start[1]
  gwas.variants$all.studies[i]  <- paste(unique(gwas.all.hits.rows$STUDY), collapse = ";    ")
  gwas.variants$risk.allele.A[i]  <- get_risk_allele_string(gwas.all.hits.rows, "A")
  gwas.variants$risk.allele.C[i]  <- get_risk_allele_string(gwas.all.hits.rows, "C")
  gwas.variants$risk.allele.G[i]  <- get_risk_allele_string(gwas.all.hits.rows, "G")
  gwas.variants$risk.allele.T[i]  <- get_risk_allele_string(gwas.all.hits.rows, "T")
  gwas.variants$risk.allele.NA[i]  <- get_risk_allele_string(gwas.all.hits.rows, "?")
  gwas.variants$pubmed.links[i]  <- paste(gwas.all.hits.rows$LINK, collapse = " ; ")
}

##################
# Getting LD Data
##################

genomes.population <- "CEU" # Using Utah 1000 genomes data

# Variants that are GWAS hits significant eQTL
eQTL.gwas.combined.rsid <- intersect(gwas.variants$RSID, eQTL.combined$SNPid)
eQTL.gwas.combined.LD.r2 <-matrix(nrow = length(eQTL.gwas.combined.rsid), ncol= length(eQTL.gwas.combined.rsid))
eQTL.gwas.combined.LD.dprime <-matrix(nrow = length(eQTL.gwas.combined.rsid), ncol= length(eQTL.gwas.combined.rsid))
colnames(eQTL.gwas.combined.LD.r2) <- colnames(eQTL.gwas.combined.LD.dprime) <- eQTL.gwas.combined.rsid
rownames(eQTL.gwas.combined.LD.r2) <- rownames(eQTL.gwas.combined.LD.dprime) <- eQTL.gwas.combined.rsid
for(i in 1:(length(eQTL.gwas.combined.rsid) - 1)){
  for(j in (i+1):length(eQTL.gwas.combined.rsid)){
    ensembl.ld.json <- read_json(paste("http://grch37.rest.ensembl.org/ld/human/pairwise/",
                                       eQTL.gwas.combined.rsid[i],"/",eQTL.gwas.combined.rsid[j],
                                       "?content-type=application/json;population_name=1000GENOMES:",
                                       "phase_3:",genomes.population,sep = ""))
    if(length(ensembl.ld.json) > 0){
      eQTL.gwas.combined.LD.r2[i,j] <- eQTL.gwas.combined.LD.r2[j, i] <- ensembl.ld.json[[1]]$r2
      eQTL.gwas.combined.LD.dprime[i,j] <- eQTL.gwas.combined.LD.dprime[j, i] <- ensembl.ld.json[[1]]$d_prime
    }
  }
}

# LD data within 500 KB of the center of the gene of interest
gene.of.interest.half <- (gene.of.interest.start + gene.of.interest.end) %/% 2
start.500Kb <- gene.of.interest.half - 250000
end.500Kb <- gene.of.interest.half + 250000
LD.info.500Kb <- read_json(paste("http://grch37.rest.ensembl.org/ld/human/region/",gene.of.interest.chrom,":",
                                 start.500Kb,"..",end.500Kb -1,"/1000GENOMES:phase_3:",genomes.population,
                                 "?content-type=application/json",sep = ""))
LD.info.500Kb.length <- length(LD.info.500Kb)
LD.info.500Kb.all.variants <- vector(length = 1e06)
for(n in 1:LD.info.500Kb.length){
  LD.info.500Kb.all.variants[n] <- (LD.info.500Kb[[n]]$variation1)
  LD.info.500Kb.all.variants[LD.info.500Kb.length + n] <- (LD.info.500Kb[[n]]$variation2)
}

LD.info.500Kb.unique.variants <- data.frame(rsid = unique(LD.info.500Kb.all.variants))
LD.info.500Kb.unique.variants$position <- numeric(dim(LD.info.500Kb.unique.variants)[1])

# Fetch Position from Ensembl using RSID
server <- "http://grch37.rest.ensembl.org"
ext <- "/variation/homo_sapiens"
i <- 1
while(i < dim(LD.info.500Kb.unique.variants)[1]){
  j <- i + 190 # Ensembl takes at most 200 requests at a time.
  if(j > dim(LD.info.500Kb.unique.variants)[1]){
    j = dim(LD.info.500Kb.unique.variants)[1]
  }
  rest.api.response <- r <- POST(paste(server, ext, sep = ""), content_type("application/json"), accept("application/json"),
                            body = paste('{ "ids" : [', paste0(LD.info.500Kb.unique.variants$rsid[i:j],
                                                               collapse = "\",\""), ' ] }', sep = "\""))
  rest.api.info <- fromJSON(toJSON(content(rest.api.response)))
  for(k in 1:length(rest.api.info)){
    LD.info.500Kb.unique.variants$position[k + i - 1] <- rest.api.info[[k]]$mappings$start[[1]]
  }
  i <- j + 1
}

# Get rid of 0s if applicable and sort ascending by position
LD.info.500Kb.unique.variants <- LD.info.500Kb.unique.variants[which(LD.info.500Kb.unique.variants$position > 0),]
LD.info.500Kb.unique.variants <- LD.info.500Kb.unique.variants[order(LD.info.500Kb.unique.variants$position),]


LD.r2.500Kb <- matrix(0,nrow = dim(LD.info.500Kb.unique.variants)[1], ncol= dim(LD.info.500Kb.unique.variants)[1])
LD.dprime.500Kb <- matrix(0,nrow = dim(LD.info.500Kb.unique.variants)[1], ncol= dim(LD.info.500Kb.unique.variants)[1])

LD.500Kb = data.frame(from = character(length(LD.info.500Kb)),
                to = character(length(LD.info.500Kb)),
                r2 = numeric(length(LD.info.500Kb)),
                dprime = numeric(length(LD.info.500Kb)),
                position.from = numeric(length(LD.info.500Kb)),
                position.to = numeric(length(LD.info.500Kb)),
                stringsAsFactors = FALSE)

rownames(LD.r2.500Kb) <-colnames(LD.r2.500Kb) <- rownames(LD.dprime.500Kb) <- colnames(LD.dprime.500Kb) <-
  LD.info.500Kb.unique.variants$rsid

j<-1
for (i in LD.info.500Kb) {
  LD.r2.500Kb[i$variation1,i$variation2] <- LD.r2.500Kb[i$variation2,i$variation1] <- as.numeric(i$r2)
  LD.dprime.500Kb[i$variation1,i$variation2] <- LD.dprime.500Kb[i$variation2,i$variation1] <-
    as.numeric(i$d_prime)
  LD.500Kb$from[j] <- i$variation1
  LD.500Kb$position.from[j] <- LD.info.500Kb.unique.variants$position[
    which(LD.info.500Kb.unique.variants$rsid == i$variation1)]
  LD.500Kb$to[j]<-i$variation2
  LD.500Kb$position.to[j] <- LD.info.500Kb.unique.variants$position[
    which(LD.info.500Kb.unique.variants$rsid == i$variation2)]
  LD.500Kb$r2[j] <- as.numeric(i$r2)
  LD.500Kb$dprime[j] <- as.numeric(i$d_prime)
  j <- j + 1
}

i<-data.frame(rsid = names(LD.r2.500Kb[gwas.variants$RSID[1],]),r2 = LD.r2.500Kb[gwas.variants$RSID[1],],
              row.names = NULL)
i <- i[which(i$r2 > 0),]
i <- i[order(i$r2, decreasing = TRUE),]

j<-data.frame(rsid = names(LD.r2.500Kb[gwas.variants$RSID[2],]),r2 = LD.r2.500Kb[gwas.variants$RSID[2],],
              row.names = NULL)
j <- j[which(j$r2 > 0),]
j<- j[order(j$r2, decreasing = TRUE),]


##############
# Gene Plot
##############

minbase <- min( as.numeric(eQTL.combined$position))
maxbase <- max( as.numeric(eQTL.combined$position))
gene.image <- makeGene(id = gene.of.interest, type = "hgnc_symbol", biomart = mart)
genesplus <- makeGeneRegion(start = minbase, end = maxbase,
                            strand = "+", chromosome = gene.of.interest.chrom, biomart=mart)
genesmin <- makeGeneRegion(start = minbase, end = maxbase,
                           strand = "-", chromosome = gene.of.interest.chrom, biomart=mart)
expres.glom <- makeSegmentation(value = as.numeric(-log10(as.numeric(eQTL.combined$pvalue[
  which(eQTL.combined$compartment == "Glom")]))), start = as.numeric(
    eQTL.combined$position[which(eQTL.combined$compartment == "Glom")]),end = as.numeric(
    eQTL.combined$position[which(eQTL.combined$compartment == "Glom")]),dp = DisplayPars(color="darkblue", lwd = 3, lty=1))
expres.tub <- makeGenericArray(intensity = as.matrix(-log10(as.numeric(eQTL.combined$pvalue[
  which(eQTL.combined$compartment == "Tub")]))), probeStart = as.numeric(
  eQTL.combined$position[which(eQTL.combined$compartment == "Tub")]), dp = DisplayPars(type = "dot", lwd = 3, pch = 5),
  trackOverlay = expres.glom)
genomeAxis <- makeGenomeAxis(add53 = TRUE, add35=TRUE)
gene.region.overlay <- makeRectangleOverlay(start = gene.of.interest.start, end = gene.of.interest.end,
                                            dp = DisplayPars(fill = "yellow", alpha = 0.2, lty = "dashed"))
legend = makeLegend(text = c('Tub','Glom', "GWAS", gene.of.interest),
                    fill = c('darkred','darkblue', "darkgreen", "lightyellow"), cex = 1)

overlays <- vector(mode="list",length = dim(gwas.variants)[1] + 1)
overlays[dim(gwas.variants)[1] + 1]<- makeRectangleOverlay(
  start = gene.of.interest.start,end = gene.of.interest.end,
  dp = DisplayPars(fill = "yellow", alpha = 0.2, lty = "dotted"), region = c(2,3))
for(i in 1:dim(gwas.variants)[1]){
  overlays[i]<- makeTextOverlay("o", xpos = as.numeric(gwas.variants$position[i]), ypos = .13,
                                     coords = "genomic", dp = DisplayPars(color = "darkgreen", cex = 1.5))
}
gdPlot(list(legend,"-log(P value)" = expres.tub, BP = genomeAxis), overlays = overlays,
       minBase = minbase, maxBase =maxbase, labelCex = 2)

# Graph "zoomed in" to 100,000 range around gene
zoom.minbase <- gene.of.interest.start - 15000
zoom.maxbase <- gene.of.interest.end + 15000
zoom.overlays <- vector(mode="list",length = dim(gwas.variants)[1] + 1)
zoom.overlays[dim(gwas.variants)[1] + 1]<- makeRectangleOverlay(
  start = gene.of.interest.start,end = gene.of.interest.end,
  dp = DisplayPars(fill = "yellow", alpha = 0.2, lty = "dotted"), region = c(2,3))
for(i in 1:dim(gwas.variants)[1]){
  zoom.overlays[i]<- makeTextOverlay("o", xpos = as.numeric(gwas.variants$position[i]), ypos = .1,
                                     coords = "genomic", dp = DisplayPars(color = "darkgreen", cex = 1.5))
}
zoom.legend = makeLegend(text = c('Tub','Glom', gene.of.interest, "GWAS"),
                    fill = c('darkred','darkblue', "lightyellow", "darkgreen"), cex = 1)
gdPlot(list(zoom.legend, "-log(P value)" = expres.tub, genomeAxis),
       overlays = zoom.overlays,
       minBase = zoom.minbase, maxBase = zoom.maxbase, labelCex = 2)

# gdPlot(list(genomeAxis, gene.image, expres.gwas), minBase = min(as.numeric(gwas.variants$position) ),
#        maxBase = as.numeric(max(gwas.variants$position) ))

get_breaks<-function(minbase, maxbase){
  lower <- floor(minbase/1000)*1000
  upper <- ceiling(maxbase/1000)*1000
  seglen <- (upper - lower) / 5
  return(seq(from = lower, to = upper, by= seglen))
}

# Creates scatterplot of eQTL hits
ggplot2::ggplot(eQTL.combined, ggplot2::aes(x = as.integer(position), y = -log10(as.numeric(eQTL.combined$pvalue)),
                                            colour = compartment)) +
  ggplot2::geom_point() +
  ggplot2:::scale_x_continuous(paste('Position on Chromosome ',gene.of.interest.chrom," (hg19)",sep=""),
                               breaks=get_breaks(minbase, maxbase)) +
  ggplot2::ylab("-log 10( P value)") +
  ggplot2::ggtitle(paste("eQTL Distribution for ", gene.of.interest, sep = ""))+
  ggplot2::theme(plot.title = ggplot2::element_text(family = "Trebuchet MS", color="#666666", face="bold", size=18, hjust=0)) +
  ggplot2::theme(axis.title = ggplot2::element_text(family = "Trebuchet MS", color="#666666", face="bold", size=14))




# Getting Start and End information for gene of interest
# ensembl.GRCh37 = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
# ensembl.GRCh38 = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# gene.of.interest.start.GRCh37 <- getBM("start_position", filters="hgnc_symbol",
#                                        values=gene.of.interest, mart=ensembl.GRCh37)[[1]]
# gene.of.interest.end.GRCh37 <- getBM("end_position", filters="hgnc_symbol",
#                                      values=gene.of.interest, mart=ensembl.GRCh37)[[1]]
# gene.of.interest.start.GRCh38 <- getBM("start_position", filters="hgnc_symbol",
#                                        values=gene.of.interest, mart=ensembl.GRCh38)[[1]]
# gene.of.interest.end.GRCh38 <- getBM("end_position", filters="hgnc_symbol",
#                                      values=gene.of.interest, mart=ensembl.GRCh38)[[1]]

# gwas.all.hits$distance.GRCh37 <- numeric(dim(gwas_variants)[1])
# gwas.all.hits$distance.GRCh38 <- numeric(dim(gwas_variants)[1])
# get_distance_from_gene <- function(n){
#   if(gwas_variants$Upstream[n]){
#     gwas_variants$distance.GRCh37[n] <<- gwas_variants$position.GRCh37[n] -
#       gene.of.interest.end.GRCh37
#     gwas_variants$distance.GRCh38[n] <<- gwas_variants$position.GRCh38[n] -
#       gene.of.interest.end.GRCh38
#   }
#   else if(gwas_variants$Downstream[n]){
#     gwas_variants$distance.GRCh37[n] <<- gene.of.interest.start.GRCh37 - gwas_variants$position.GRCh37[n]
#     gwas_variants$distance.GRCh38[n] <<- gene.of.interest.start.GRCh38 - gwas_variants$position.GRCh38[n]
#   }
#   else{
#     return(NULL)
#   }
# }
# sapply(1:dim(gwas_variants)[1], get_distance_from_gene)
# gwas_variants$Closest <- sapply(1:dim(gwas_variants)[1], function(n){
#   return(gwas_variants$genomicContexts[[n]]$isClosestGene[which(gwas_variants$genomicContexts[[n]]
#                                                                $gene$geneName == gene.of.interest)[1]][[1]])})
# gwas_variants$closest.genes <- sapply(1:dim(gwas_variants)[1], function(n){
#   return(gsub(",",";", toString(unique(gwas_variants$genomicContexts[[n]]$gene
#                                        $geneName[which(gwas_variants$genomicContexts[[n]]
#                                                        $isClosestGene == TRUE)]))))
# })

zoom.ld <- LD.500Kb[which(LD.500Kb$position.from < zoom.maxbase & LD.500Kb$position.from > zoom.minbase &
                            LD.500Kb$position.to < zoom.maxbase & LD.500Kb$position.to > zoom.minbase ),]


ld.eqtl.overlap <- LD.500Kb[which(LD.500Kb$from %in% eQTL.combined$SNPid & LD.500Kb$to %in% eQTL.combined$SNPid),]

zoom.ld.eqtl.overlap <- ld.eqtl.overlap[which(ld.eqtl.overlap$position.from < zoom.maxbase &
                                                ld.eqtl.overlap$position.from > zoom.minbase &
                                                ld.eqtl.overlap$position.to < zoom.maxbase &
                                                ld.eqtl.overlap$position.to > zoom.minbase ),]

zoom.ld.eqtl.overlap.r2 <- matrix(0, nrow = length(unique(
  union(zoom.ld.eqtl.overlap$to, zoom.ld.eqtl.overlap$from))),
  ncol = length(unique(union(zoom.ld.eqtl.overlap$to, zoom.ld.eqtl.overlap$from))),
  dimnames = list(unique(union(zoom.ld.eqtl.overlap$to, zoom.ld.eqtl.overlap$from))[order(LD.info.500Kb.unique.variants$position[
    match(unique(union(zoom.ld.eqtl.overlap$to, zoom.ld.eqtl.overlap$from)), LD.info.500Kb.unique.variants$rsid)])],
    unique(union(zoom.ld.eqtl.overlap$to, zoom.ld.eqtl.overlap$from))[order(LD.info.500Kb.unique.variants$position[
      match(unique(union(zoom.ld.eqtl.overlap$to, zoom.ld.eqtl.overlap$from)), LD.info.500Kb.unique.variants$rsid)])]))

for(i in 1:dim(zoom.ld.eqtl.overlap)[1]){
  zoom.ld.eqtl.overlap.r2[zoom.ld.eqtl.overlap$from[i], zoom.ld.eqtl.overlap$to[i]]<-
    zoom.ld.eqtl.overlap.r2[zoom.ld.eqtl.overlap$to[i],zoom.ld.eqtl.overlap$from[i]] <-
                                zoom.ld.eqtl.overlap$r2[i]
}

annot <- HeatmapAnnotation("-log10(Pvalue)" = anno_points(eqtl.combined.tub$Mlog[
  match(rownames(zoom.ld.eqtl.overlap.r2),eqtl.combined.tub$SNPid)],axis = TRUE, axis_side = "right",
  which = "column",annotation_height = unit(c(5), "cm"), ylim = c(0, 5)), show_annotation_name = TRUE,
  annotation_name_rot = 0, annotation_name_offset = unit(8, "mm"))

zoom.ld.eqtl.overlap.r2[lower.tri(zoom.ld.eqtl.overlap.r2)] = 0
col_fun = colorRamp2(c( 0, 1), c("white", "darkred"), transparency = 0.5)
Heatmap(zoom.ld.eqtl.overlap.r2, name = "LD R2", col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = TRUE, show_column_dend = TRUE, top_annotation = annot,
        top_annotation_height = unit(2, "cm"))


chordDiagram(zoom.ld.eqtl.overlap.r2, col= col_fun(zoom.ld.eqtl.overlap.r2),
             annotationTrack = c( "grid"),annotationTrackHeight = c(0.03, 0.01),
             grid.col = NA, grid.border = "black",
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(zoom.ld.eqtl.overlap.r2))))))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))}, bg.border = NA)

# eqtl.combined.tub <- eQTL.combined[which(eQTL.combined$compartment == "Tub"),]
# eqtl.combined.tub$pvalue <- as.numeric(as.character(eqtl.combined.tub$pvalue))
# eqtl.combined.tub$Mlog <- -log10(eqtl.combined.tub$pvalue)
# eqtl.combined.tub.length <- dim(eqtl.combined.tub)[1]
# eqtl.combined.glom <- eQTL.combined[which(eQTL.combined$compartment == "Glom"),]
# eqtl.combined.glom$pvalue <- as.numeric(as.character(eqtl.combined.glom$pvalue))
# eqtl.combined.glom$Mlog <- -log10(eqtl.combined.glom$pvalue)
# eqtl.combined.glom.length <- dim(eqtl.combined.glom)[1]
# bed.eqtl.tub <- data.frame(chr = character(eqtl.combined.tub.length),
#                             start = as.integer(eqtl.combined.tub$position),
#                             end = as.integer(eqtl.combined.tub$position),
#                             value = eqtl.combined.tub$pvalue, stringsAsFactors = FALSE)
# bed.eqtl.glom <- data.frame(chr = character(eqtl.combined.glom.length),
#                            start = as.integer(eqtl.combined.glom$position),
#                            end = as.integer(eqtl.combined.glom$position),
#                            value = eqtl.combined.glom$pvalue, stringsAsFactors = FALSE)
# bed.eqtl.glom$chr <- bed.eqtl.tub$chr <- paste("chr", gene.of.interest.chrom, sep = "")
#
# bed.list.eqtl <- list(bed.eqtl.glom, bed.eqtl.tub)
#
circos.clear()
basetrack = data.frame(
  name  = paste("Chrom. ",gene.of.interest.chrom, sep = ""),
  start = c(minbase),
  end   = c(maxbase))
circos.par("gap.degree" = rep(5),"start.degree" = 90)
circos.initializeWithIdeogram(chromosome.index = "chr1")
circos.genomicInitialize(basetrack)
circos.genomicTrackPlotRegion(bed.list.eqtl)

zoom.ld.eqtl.overlap.positions <- LD.info.500Kb.unique.variants$position[
  match(unique(union(zoom.ld.eqtl.overlap$to, zoom.ld.eqtl.overlap$from)),
        LD.info.500Kb.unique.variants$rsid)]
edgelist <- cbind(zoom.ld.eqtl.overlap$from, zoom.ld.eqtl.overlap$to)
arcplot(edgelist, vertices = unique(union(zoom.ld.eqtl.overlap$to, zoom.ld.eqtl.overlap$from)),
        ordering = unique(union(zoom.ld.eqtl.overlap$to, zoom.ld.eqtl.overlap$from))[
          order(zoom.ld.eqtl.overlap.positions)])

