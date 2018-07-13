###################
# Kate Sanders
# BSRP 2018
###################


library(gwascat)
library(readxl)


# Run this only if this value doesn't already exist. It takes a while.
# current_gwascat <- makeCurrentGwascat(genome = "GRCh37")

# Enter the filepath you would like all documents produced to go to
gwas.filepath <- "/Users/ksanders/Documents/"

# Put in path to Table S3 from http://science.sciencemag.org/content/360/6390/758/tab-figures-data
single.cell.all.genes <- read_excel("/Users/ksanders/Downloads/aar2131_Table_S3.xlsx",
                                    sheet= "summary_for_all_genes")[-c(1), ]

# Renaming columns
colnames(single.cell.all.genes)[colnames(single.cell.all.genes)=="X__1"] <- "gene"

sapply(2:dim(single.cell.all.genes)[2], function(n){
  name.paren <- paste(strsplit(names(single.cell.all.genes)[n], split = " ")[[1]][-1], collapse = "-")
  colnames(single.cell.all.genes)[colnames(single.cell.all.genes)==names(single.cell.all.genes)[n]] <<-
    substring(name.paren, first = 2, last = str_length(name.paren) - 1)
})

s<- c("Chronic kidney disease","Chronic kidney disease (chronic kidney disease vs normal or mildly reduced eGFR) in type 1 diabetes","Chronic kidney disease (severe chronic kidney disease vs normal kidney function) in type 1 diabetes","Chronic kidney disease and serum creatinine levels","Diabetic kidney disease","Fractional excretion of metabolites in chronic kidney disease","Glomerular filtration rate in chronic kidney disease","Immunoglobulin light chain (AL) amyloidosis (heart and kidney involvement)","Immunoglobulin light chain (AL) amyloidosis (kidney involvement)","Kidney cancer","Kidney disease (early and late stages) in type 1 diabetes","Kidney disease (early stage) in type 1 diabetes","Kidney disease (end stage renal disease vs non-end stage renal disease) in type 1 diabetes","Kidney disease (end stage renal disease vs normoalbuminuria) in type 1 diabetes","Kidney disease (late stage) in type 1 diabetes","Kidney function decline traits","Kidney stones","Proteinuria and chronic kidney disease","Proteinuria in chronic kidney disease","Renal function and chronic kidney disease","Serum metabolite concentrations in chronic kidney disease","Serum metabolite ratios in chronic kidney disease","Tacrolimus trough concentration in kidney transplant patients","Urinary metabolite concentrations in chronic kidney disease","Urinary metabolite ratios in chronic kidney disease")
