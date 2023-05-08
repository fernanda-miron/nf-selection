## Using rehh
## Upload package
library("rehh")
library("vcfR")
library("ggplot2")
library("qqman")
library("dplyr")

## Read args from command line
args <- commandArgs(trailingOnly = T)

# Uncomment for debbuging
# Comment for production mode only
#args[1] <- "./final.annotated_22.phased.with.ref.vcf"
#args[2] <- "ihs_chromosome.tsv"

## Place args into named object
file_ihs <- args[1]
export_ihs <- args[2]

## First step
haplohh <- data2haplohh(hap_file = file_ihs,remove_multiple_markers = T,
                        vcf_reader = "vcfR")

## Second step
ihh <- scan_hh(haplohh)

## Third step
ihs <- ihh2ihs(ihh)

## Extract ihs
final_ihs <- ihs$ihs

## Removing nas and making things numeric
final2 <- final_ihs %>% 
  mutate(CHR = as.numeric(final_ihs$CHR))

## remove nas
final3 <- na.omit(final2)

## creating column of pvalues
final4 <- final3 %>% 
  mutate(pvals = 10^-LOGPVALUE)

##
## export 
write.table(final4,file = export_ihs,sep = "\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
