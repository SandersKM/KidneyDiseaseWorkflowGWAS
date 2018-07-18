###################
# Kate Sanders
# BSRP 2018
###################

# source("http://bioconductor.org/biocLite.R")
# biocLite("GenomeGraphs")
library(jsonlite)
library(rentrez)
library(gwascat)
library(magrittr)
library(GenomeGraphs)
library(ensembldb)
library(ggplot2)
library(LDheatmap)

# Run this only if this value doesn't already exist. It takes a while.
# current_gwascat <- makeCurrentGwascat(genome = "GRCh37")

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
gene.of.interest.info <- getBM(c("start_position", "end_position", "strand"), filters="hgnc_symbol",
                               values=gene.of.interest, mart=mart)
gene.ofinterest.start <- gene.of.interest.info$start_position
gene.ofinterest.strand <- gene.of.interest.info$strand
gene.ofinterest.end <- gene.of.interest.info$end_position

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
    r$STRONGEST.SNP.RISK.ALLELE, first = str_length(r$STRONGEST.SNP.RISK.ALLELE[1])) == a),]
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

# Variants that are GWAS hits significant eQTL
eQTL.gwas.combined.rsid <- intersect(gwas.variants$RSID, eQTL.combined$SNPid)
eQTL.gwas.combined.LD.r2 <-matrix(nrow = length(eQTL.gwas.combined.rsid), ncol= length(eQTL.gwas.combined.rsid))
eQTL.gwas.combined.LD.dprime <-matrix(nrow = length(eQTL.gwas.combined.rsid), ncol= length(eQTL.gwas.combined.rsid))
colnames(eQTL.gwas.combined.LD.r2) <- colnames(eQTL.gwas.combined.LD.dprime) <- eQTL.gwas.combined.rsid
rownames(eQTL.gwas.combined.LD.r2) <- rownames(eQTL.gwas.combined.LD.dprime) <- eQTL.gwas.combined.rsid
genomes.population <- "CEU" # Using Utah 1000 genomes data
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

##############
# Gene Plot
##############

minbase <- min( as.numeric(eQTL.combined$position))
maxbase <- max( as.numeric(eQTL.combined$position))
gene.image <- makeGene(id = gene.of.interest, type = "hgnc_symbol", biomart = mart)
gene.of.interest.chrom <- eQTL.combined$chrom[1]
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
gene.region.overlay <- makeRectangleOverlay(start = gene.ofinterest.start, end = gene.ofinterest.end,
                                            dp = DisplayPars(fill = "yellow", alpha = 0.2, lty = "dashed"))
legend = makeLegend(text = c('Tub','Glom', gene.of.interest),
                    fill = c('darkred','darkblue', "lightyellow"), cex = 1)
gdPlot(list(genesplus,genomeAxis,genesmin, "-log(P value)" = expres.tub, legend), overlays = gene.region.overlay,
       minBase = minbase, maxBase =maxbase, labelCex = 2)

# Graph "zoomed in" to 100,000 range around gene
zoom.minbase <- gene.of.interest.start - 100000
zoom.maxbase <- gene.of.interest.end + 100000
zoom.overlays <- vector(mode="list",length = dim(gwas.variants)[1] + 1)
zoom.overlays[dim(gwas.variants)[1] + 1]<- makeRectangleOverlay(
  start = gene.ofinterest.start,end = gene.ofinterest.end,
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


