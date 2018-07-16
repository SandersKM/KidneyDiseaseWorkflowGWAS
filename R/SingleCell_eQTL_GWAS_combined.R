###################
# Kate Sanders
# BSRP 2018
###################

library(pheatmap)
library(gwascat)
library(scater)
library(readxl)
library(stringi)

# Replace the following with the specific gwas traits of interest
gwas.traits.of.interest <- c("Chronic kidney disease (chronic kidney disease vs normal
                             or mildly reduced eGFR) in type 1 diabetes","Chronic kidney disease
                             (severe chronic kidney disease vs normal kidney function) in type 1 diabetes")

# Run this only if this value doesn't already exist. It takes a while.
# current_gwascat <- makeCurrentGwascat(genome = "GRCh37")

# Enter the filepath you would like all documents produced to go to
gwas.filepath <- "/Users/ksanders/Documents/"

# Put in path to Table S3 from http://science.sciencemag.org/content/360/6390/758/tab-figures-data
single.cell.all.genes <- read_excel("/Users/ksanders/Downloads/aar2131_Table_S3.xlsx",
                                    sheet= "summary_for_all_genes")[-c(1), ]

# Renaming columns
colnames(single.cell.all.genes)[colnames(single.cell.all.genes)=="X__1"] <- "mouse.gene"

sapply(2:dim(single.cell.all.genes)[2], function(n){
  name.paren <- paste(strsplit(names(single.cell.all.genes)[n], split = " ")[[1]][-1], collapse = "-")
  colnames(single.cell.all.genes)[colnames(single.cell.all.genes)==names(single.cell.all.genes)[n]] <<-
    substring(name.paren, first = 2, last = stri_length(name.paren) - 1)
})

# Converting Mouse genes to Human Genes
mart <- get0("mart", ifnotfound=useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37))
mouse.mart <- get0("mouse.mart", ifnotfound=useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", GRCh=37))
single.cell.all.genes <- cbind(human.gene = "", single.cell.all.genes)
mouse.to.human <- function(x){
  genesV2 <- getLDS(attributes = "mgi_symbol", filters = "mgi_symbol", values = x , mart = mouse.mart,
                    attributesL = "hgnc_symbol", martL = mart, uniqueRows=T)
  return(unique(genesV2[, 2]))
}
single.cell.all.genes$human.gene <- sapply(single.cell.all.genes$mouse.gene, mouse.to.human)
# search for gwas hits
gwas.hits.by.trait <- subsetByTraits(current_gwascat, tr = gwas.traits.of.interest)

