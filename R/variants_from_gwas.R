###################
# Kate Sanders
# BSRP 2018
###################

# library("car")
# source("http://bioconductor.org/biocLite.R")
# biocLite("GenomeGraphs")
library(httr)
library(xml2)
library(jsonlite)
library(rentrez)
library(gwascat)
library(magrittr)
library(GenomeGraphs)
library(ensembldb)
library(eqtl)

# Run this only if this value doesn't already exist. It takes a while.
current_gwascat <- makeCurrentGwascat(genome = "GRCh37")

# Enter the HGNC symbol of the gene of interest below
gene.of.interest <- "MUC1"

# Enter the filepath you would like all documents produced to go to
gwas.filepath <- "/Users/ksanders/Documents/"

# This will bring you to the nephvs eQTL webpage
browseURL(paste("http://eqtl.nephvs.org/searchResult/", gene.of.interest, sep = ""))
# Enter the path to the location of these downloaded files
filepath <- "/Users/ksanders/Desktop/"
nephQTL.glom <- read.csv(paste(filepath, "glom_MatrixEQTL_", gene.of.interest,".csv", sep = ""),
                         header = TRUE, sep = ",")
nephQTL.tub <- read.csv(paste(filepath, "tub_MatrixEQTL_", gene.of.interest,".csv", sep = ""),
                        header = TRUE, sep = ",")

# Copy and paste the tubule eQTL for gene of interest from http://18.217.22.69/eqtl
Susztak.tub <- read.csv(paste(filepath,"Susztak.eQTL.Tubule.MUC1",".csv", sep = ""),
                        header = TRUE, sep = ",")

# Make general table of eQTL positions/values
total.rows <- dim(nephQTL.glom)[1] + dim(nephQTL.tub)[1]
eQTL.combined <- data.frame(SNPid = character(total.rows), chrom = character(total.rows),
                            position = character(total.rows), ref = character(total.rows),
                            alt = character(total.rows), pvalue = character(total.rows),
                            beta = character(total.rows), compartment = character(total.rows),
                            source = character(total.rows), stringsAsFactors = FALSE)

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
# sapply(1:dim(Susztak.tub)[1], rowstart = dim(nephQTL.glom)[1] + dim(nephQTL.tub)[1], function(n, rowstart){
#   i <- rowstart + n
#   eQTL.combined$SNPid[i] <<- toString(Susztak.tub$SNP[n])
#   eQTL.combined$chrom[i]<<-toString(Susztak.tub$Chr[n])
#   eQTL.combined$position[i]<<-toString(Susztak.tub$Loc[n])
#   eQTL.combined$ref[i]<<-toString(Susztak.tub$Ref.Allele[n])
#   eQTL.combined$alt[i]<<-toString(Susztak.tub$Alt.Allele[n])
#   eQTL.combined$pvalue[i]<<-toString(Susztak.tub$Pval[n])
#   eQTL.combined$beta[i]<<-toString(Susztak.tub$Beta[n])
#   eQTL.combined$compartment[i]<<- "Tub"
#   eQTL.combined$source[i]<<- "Susztak"
# })

# Get combined sheet with GWAS and eQTL



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

reported_gene_rows <- as.integer(get_reported_gene_of_interest(length(
  current_gwascat$`REPORTED GENE(S)`))[[1]])
mapped_gene_rows <- which(current_gwascat$MAPPED_GENE == "MUC1")
combined_gene_rows <- unique(append(mapped_gene_rows, reported_gene_rows))
# Data frame with combined gene rows with gwas hits
gwas.specific <- as.data.frame(current_gwascat[combined_gene_rows])

unique.strongest.risk.allele <- unique(gwas.specific$STRONGEST.SNP.RISK.ALLELE)
gwas.txt.file <- gene.of.interest
for(i in unique.strongest.risk.allele){
  allele.rows <- which(gwas.specific$STRONGEST.SNP.RISK.ALLELE == i)
  gwas.txt.file <- paste(gwas.txt.file, paste(
    "RSid: ", i, "; Position (hg19): ",gwas.specific$seqname[allele.rows[1]], ":",
    gwas.specific$start[allele.rows[1]], sep = ""), sep = "\n")
  n <- 1
  for(j in allele.rows){
    row <- gwas.specific[j,]
    entrez.var <- entrez_summary(db="pubmed", id=row$PUBMEDID)
    gwas.txt.file <- paste(gwas.txt.file, "\n",
      paste("\t", n, ") Study: ", entrez.var$title, sep = ""),
      paste("First Author: ", entrez.var$firstauthor, ";   Last Author: ", entrez.var$lastauthor, sep = ""),
      paste("Journal: ", entrez.var$fulljournalname,";    Date Published: ", entrez.var$pubdate, sep=""),
      paste(entrez.var$articleids$idtype, ":", entrez.var$articleids$value,collapse = " ; "),
      paste("Disease Trait: ", row$DISEASE.TRAIT, ";    Mapped Trait: ", row$MAPPED.TRAIT),
      paste("Reported / Mapped Genes: ", gsub(",","; ",paste(unique(c(
        row$REPORTED.GENE.S., row$MAPPED_GENE)), collapse = "; " )), sep = ""),
      paste("Initial Sample Size: ", row$INITIAL.SAMPLE.SIZE,
            ";     Replication Sample Size: ", row$REPLICATION.SAMPLE.SIZE, sep = ""),
      paste("Risk Allele Frequency: ", row$RISK.ALLELE.FREQUENCY, sep = ""),
      paste("Platform: ", strsplit(row$PLATFORM..SNPS.PASSING.QC., split = " ")[[1]][1],
           ";    SNPs Passing WC: ", paste(strsplit(row$PLATFORM..SNPS.PASSING.QC., split = " ")[[1]][-1],
                                           collapse = ""), sep= ""),
      paste("P-Value: ", row$P.VALUE, ";    -log_10(P-Value): ", row$PVALUE_MLOG, ";   Odds Ratio / BETA: ",
            row$OR.or.BETA, sep=""),
      sep = "\n\t")
    n <- n+1
  }
}
gwas.txt.file.write <- file(paste(gwas.filepath, "gwas_", gene.of.interest, ".txt",sep = ""))
writeLines(gwas.txt.file, gwas.txt.file.write)
close(gwas.txt.file.write)

gwas.variants <- data.frame(RSID = unique(gwas.specific$SNPS))

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
  allele.rows <- gwas.specific.rows[which(substring(
    r$STRONGEST.SNP.RISK.ALLELE, first = str_length(r$STRONGEST.SNP.RISK.ALLELE[1])) == a),]
  return(paste("Disease Trait: ", allele.rows$DISEASE.TRAIT,"; P-Value: ", allele.rows$P.VALUE,
               "; -log_10(P-Value): ", allele.rows$PVALUE_MLOG, "; Odds Ratio / BETA: ", allele.rows$OR.or.BETA,
               "; 95% CI: ", allele.rows$X95..CI..TEXT., collapse = " | "))
}
for(i in 1:dim(gwas.variants)[1]){
  gwas.specific.rows <- gwas.specific[which(gwas.specific$SNPS == gwas.variants$RSID[i]),]
  gwas.variants$chrom[i] <- gwas.specific.rows$CHR_ID[1]
  gwas.variants$position[i] <- gwas.specific.rows$start[1]
  gwas.variants$all.studies[i]  <- paste(unique(gwas.specific.rows$STUDY), collapse = ";    ")
  gwas.variants$risk.allele.A[i]  <- get_risk_allele_string(gwas.specific.rows, "A")
  gwas.variants$risk.allele.C[i]  <- get_risk_allele_string(gwas.specific.rows, "C")
  gwas.variants$risk.allele.G[i]  <- get_risk_allele_string(gwas.specific.rows, "G")
  gwas.variants$risk.allele.T[i]  <- get_risk_allele_string(gwas.specific.rows, "T")
  gwas.variants$risk.allele.NA[i]  <- get_risk_allele_string(gwas.specific.rows, "?")
  gwas.variants$pubmed.links[i]  <- paste(gwas.specific.rows$LINK, collapse = " ; ")
}


# Gene Plot
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# plusStrand <- makeGeneRegion(chromosome = as.numeric(eQTL.combined$chrom[1]),
#                              start =min( as.numeric(eQTL.combined$position)),
#                              end = max( as.numeric(eQTL.combined$position)), strand = "+", biomart = mart)
minStrand <- makeGeneRegion(chromosome = as.numeric(eQTL.combined$chrom[1]),
                            start =min( as.numeric(eQTL.combined$position)),
                            end = max( as.numeric(eQTL.combined$position)),
                            strand = "-", biomart = mart)
ideogram <- makeIdeogram(as.numeric(eQTL.combined$chrom[1]))
genomeAxis <- makeGenomeAxis(add53 = TRUE, add35 = TRUE)
gdPlot(list(ideogram,genomeAxis, minStrand, gene.image),
       minBase = min( as.numeric(eQTL.combined$position)),
       maxBase = max( as.numeric(eQTL.combined$position)))


gene.image <- makeGene(id = gene.of.interest, type = "hgnc_symbol", biomart = mart)
gdPlot(gene.image)

data("exampleData", package="GenomeGraphs")

gene.of.interest.chrom <- eQTL.combined$chrom[1]

genesplus <- makeGeneRegion(start = minbase, end = maxbase,
                            strand = "+", chromosome = gene.of.interest.chrom, biomart=mart)
genesmin <- makeGeneRegion(start = minbase, end = maxbase,
                           strand = "-", chromosome = gene.of.interest.chrom, biomart=mart)

ideog <- makeIdeogram(chromosome = as.numeric(gene.of.interest.chrom))
expres <- makeGenericArray(intensity = as.matrix(eQTL.combined$pvalue), probeStart = as.numeric(
  eQTL.combined$position),dp = DisplayPars(color="darkred", type="point"))
genomeAxis <- makeGenomeAxis(add53 = TRUE, add35=TRUE)

gdPlot(list(a=ideog,b = expres,d=genesplus,e=genomeAxis,f=genesmin),
       minBase = minbase, maxBase =maxbase, labelCex = 2)


# Creates
minbase <- min( as.numeric(eQTL.combined$position))
maxbase <- max( as.numeric(eQTL.combined$position))
scatter.eQTL <- ggplot2::ggplot(eQTL.combined, ggplot2::aes(x = as.integer(position), y = -log10(as.numeric(eQTL.combined$pvalue)),
                                            colour = compartment)) +
  ggplot2::geom_point() +
  ggplot2:::scale_x_continuous(paste('Position on Chromosome ',gene.of.interest.chrom," (hg19)",sep=""),
                               breaks=get_breaks(minbase, maxbase)) +
  ggplot2::ylab("-log 10( P value)") +
  ggplot2::ggtitle(paste("eQTL Distribution for ", gene.of.interest, sep = ""))+
  ggplot2::theme(plot.title = ggplot2::element_text(family = "Trebuchet MS", color="#666666", face="bold", size=18, hjust=0)) +
  ggplot2::theme(axis.title = ggplot2::element_text(family = "Trebuchet MS", color="#666666", face="bold", size=14))

get_breaks<-function(minbase, maxbase){
  lower <- floor(minbase/1000)*1000
  upper <- ceiling(maxbase/1000)*1000
  seglen <- (upper - lower) / 5
  return(seq(from = lower, to = upper, by= seglen))
}

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


