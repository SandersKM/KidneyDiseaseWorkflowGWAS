###################
# Kate Sanders
# BSRP 2018
###################

library(httr)
library(xml2)
library(jsonlite)

# Enter the HGNC symbol of the gene of interest below
gene.of.interest <- "MUC1"

# Access the GWAS cataloge records through their API
gwas.api.url <- paste("https://www.ebi.ac.uk/gwas/labs/rest/api/singleNucleotidePolymorphisms/search/findByGene?geneName=",
                      gene.of.interest, sep = "")
gwas.api.content <- fromJSON(toJSON(content(GET(gwas.api.url))))
gwas_variants <- data.frame(rsID = unlist(gwas.api.content$`_embedded`$singleNucleotidePolymorphisms$rsId))
gwas_variants$functionalClass <- unlist(gwas.api.content$`_embedded`$singleNucleotidePolymorphisms$functionalClass)
gwas_variants$genomicContexts <- sapply(1:dim(gwas_variants)[1], function(n){
  gwas.api.content$`_embedded`$singleNucleotidePolymorphisms$genomicContexts[n]})
gwas_variants$chromosome <- sapply(1:dim(gwas_variants)[1], function(n){
  gwas_variants$genomicContexts[n][[1]]$location$chromosomeName[[1]]})
gwas_variants$position.GRCh38 <- sapply(1:dim(gwas_variants)[1], function(n){
  gwas_variants$genomicContexts[n][[1]]$location$chromosomePosition[[1]]})
gwas_variants$position.GRCh37 <- numeric(dim(gwas_variants)[1])
gwas_variants$ancestral.allele <- numeric(dim(gwas_variants)[1])
gwas_variants$minor.allele <- numeric(dim(gwas_variants)[1])
# Minor Allele Frequency (MAF) refers to the lowest allele frequency of a sequence
# variant (such as a SNP). The global MAF is calculated using the allele frequences
# across all 1000 Genomes Phase I populations. - Ensembl.org
gwas_variants$MAF <- numeric(dim(gwas_variants)[1])
gwas_variants$most.severe.consequence <- character(dim(gwas_variants)[1])
gwas_variants$synonyms <- character(dim(gwas_variants)[1])
# Ensembl API imports Variants (including SNPs and indels)from dbSNP
variant_from_ensembl <- function(n){
  cont<- fromJSON(toJSON(content(GET(paste("http://grch37.rest.ensembl.org/variation/human/",
                                           gwas_variants$rsID[n],"?content-type=application/json",
                                           sep = "")))))
  gwas_variants$position.GRCh37[n] <<- cont$mappings$start
  gwas_variants$ancestral.allele[n] <<- cont$ancestral_allele
  gwas_variants$minor.allele[n] <<- cont$minor_allele
  gwas_variants$MAF[n] <<- cont$MAF
  gwas_variants$most.severe.consequence[n] <<- cont$most_severe_consequence
  gwas_variants$synonyms[n] <<- gsub(",", ";", toString(cont$synonyms))
  }
sapply(1:dim(gwas_variants)[1],variant_from_ensembl)



