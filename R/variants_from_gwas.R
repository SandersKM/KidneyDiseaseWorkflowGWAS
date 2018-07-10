###################
# Kate Sanders
# BSRP 2018
###################

library(httr)
library(xml2)
library(jsonlite)
library(biomaRt)
library(gwascat)
library(magrittr)

# Run this only if this value doesn't already exist. It takes a while.
current_gwascat <- makeCurrentGwascat(genome = "GRCh37")

# Enter the HGNC symbol of the gene of interest below
gene.of.interest <- "MUC1"

# This will bring you to the nephvs eQTL webpage
browseURL(paste("http://eqtl.nephvs.org/searchResult/", gene.of.interest, sep = ""))
# Enter the path to the location of these downloaded files
filepath <- "/Users/ksanders/Desktop/"
nephQTL.glom <- read.csv(paste(filepath, "glom_MatrixEQTL_", gene.of.interest,".csv", sep = ""),
                         header = TRUE, sep = ",")
nephQTL.tub <- read.csv(paste(filepath, "tub_MatrixEQTL_", gene.of.interest,".csv", sep = ""),
                        header = TRUE, sep = ",")


# Get GWAS results for variants reported or mapped to Gene of Interest
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

reported_gene_rows <- as.integer(get_reported_gene_of_interest(length(current_gwascat$`REPORTED GENE(S)`))[[1]])
mapped_gene_rows <- which(current_gwascat$MAPPED_GENE == "MUC1")
combined_gene_rows <- unique(append(mapped_gene_rows, reported_gene_rows))
# Data frame with combined gene rows with gwas hits
gwas.specific <- as.data.frame(current_gwascat[combined_gene_rows])




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

# gwas.specific$distance.GRCh37 <- numeric(dim(gwas_variants)[1])
# gwas.specific$distance.GRCh38 <- numeric(dim(gwas_variants)[1])
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


