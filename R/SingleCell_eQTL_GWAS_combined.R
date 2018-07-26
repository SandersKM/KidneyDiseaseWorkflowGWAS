################### Kate Sanders BSRP 2018

# source('http://bioconductor.org/biocLite.R') biocLite('biomaRt')
library(gwascat)
library(scater)
library(readxl)
library(stringi)
library(biomaRt)

# Replace the following with the specific gwas traits of interest
gwas.traits.of.interest <- c("Chronic kidney disease (chronic kidney disease vs normal or mildly reduced eGFR) in type 1 diabetes", "Chronic kidney disease (severe chronic kidney disease vs normal kidney function) in type 1 diabetes")

# Enter the filepath you would like all documents produced to go to
gwas.filepath <- "/Users/ksanders/Documents/"

# Put in path to Table S3 from http://science.sciencemag.org/content/360/6390/758/tab-figures-data
single.cell.all.genes <- read_excel("/Users/ksanders/Downloads/aar2131_Table_S3.xlsx", sheet = "summary_for_all_genes")[-c(1), ]

# Renaming columns
colnames(single.cell.all.genes)[colnames(single.cell.all.genes) == "X__1"] <- "mouse.gene"

sapply(2:dim(single.cell.all.genes)[2], function(n) {
    name.paren <- paste(strsplit(names(single.cell.all.genes)[n], split = " ")[[1]][-1], collapse = "-")
    colnames(single.cell.all.genes)[colnames(single.cell.all.genes) == names(single.cell.all.genes)[n]] <<- substring(name.paren, first = 2, 
        last = stri_length(name.paren) - 1)
    single.cell.all.genes[[n]] <<- as.numeric(single.cell.all.genes[[n]])
})

# Loading current GWAS Cateloge information.  Do this every few months, but be warned it take a long time to load
if (!exists("current_gwascat")) {
    current_gwascat <- makeCurrentGwascat(genome = "GRCh37")
}

# Getting Human genome information
if (!exists("mart")) {
    mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37)
}
# Getting Mouse genome information
if (!exists("mouse.mart")) {
    mouse.mart <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", GRCh = 37)
}
# converting mouse gene names to human gene names
if (!exists("genesV2")) {
    genesV2 <- getLDS(attributes = "mgi_symbol", filters = "mgi_symbol", values = single.cell.all.genes$mouse.gene, mart = mouse.mart, 
        attributesL = "hgnc_symbol", martL = mart, uniqueRows = T)
}
names(genesV2) <- c("mouse.gene", "human.gene")
single.cell.all.genes <- merge(x = single.cell.all.genes, y = genesV2, by = "mouse.gene", all.x = T, all.y = F)
single.cell.all.genes <- single.cell.all.genes[, c(ncol(single.cell.all.genes), 1:(ncol(single.cell.all.genes) - 1))]

# search for gwas hits
gwas.hits.by.trait <- subsetByTraits(current_gwascat, tr = gwas.traits.of.interest)

gwas.mapped.genes <- union(gwas.hits.by.trait$MAPPED_GENE, gwas.hits.by.trait$`REPORTED GENE(S)`)
gwas.mapped.genes <- unique(unlist(strsplit(gsub(",", " -", gwas.mapped.genes), split = " - ")))
gwas.mapped.genes <- gwas.mapped.genes[gwas.mapped.genes != "NR"]



