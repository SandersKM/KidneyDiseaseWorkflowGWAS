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

# TODO - Convert position to GRCh37 using http://grch37.rest.ensembl.org/documentation/info/variation_id
