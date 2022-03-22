## This script takes ihs and pbs data and makes a circus plot

library("BioCircos")
library("dplyr")
library("vroom")
library("purrr")

## Read args from command line
args <- commandArgs(trailingOnly = T)

# Uncomment for debbuging
# Comment for production mode only
#args[1] <- "pbs.tsv"
#args[2] <- "final_ihs.tsv"
#args[3] <- 0.2
#args[4] <- 2

## Place args into named object
pbs_file <- args[1]
ihs_file <- args[2]
p_cut <- args[1]
i_cut <- args[2]

## Up data frame 
pbs <- vroom(pbs_file)
ihs <- vroom(ihs_file)

## Remove inf
winf <- pbs[!is.infinite(rowSums(pbs)),]

## Filter dataframe to only keep snps with higher PBS value
significative_pbs <- winf %>% 
  filter(PBS_value > 0.8)

## Filter dataframe to only keep snps with higher iHS value
significative_ihs <- ihs %>% 
  filter(Std_iHS > 5)

## Filter to get SNPs from iHS and PBS
## First, lets get the top 1% of PBS
arrange_PBS <- winf[with(winf, order(-PBS_value)),]

## Get the 1% 
top_PBS <- head(arrange_PBS, (nrow(arrange_PBS)*0.01))

## Filter to get SNPs from iHS and PBS
## Now from iHS
arrange_iHS <- ihs[with(ihs, order(-Std_iHS)),]

## Get the 1% 
top_ihs <- head(arrange_iHS, (nrow(arrange_iHS)*0.01))
names(top_ihs)[names(top_ihs) == 'chr'] <- 'CHROM'
names(significative_ihs)[names(significative_ihs) == 'chr'] <- 'CHROM'

## Merge top PBS and top iHS
merged <- left_join(top_ihs, top_PBS, by = c("POS", "CHROM"))

## Remove NAs
final <- merged[complete.cases(merged), ]

## Write table
write.table(x = final, file = "pbs_vs_ihs.tsv", sep = "\t", 
            col.names = T, row.names = F, quote = F)

## Add column with value and label
final2 <- final %>% 
  mutate(value = Std_iHS + PBS_value) %>% 
  mutate(label = paste("PBS =", PBS_value, "iHS =", Std_iHS)) %>% 
  filter(Std_iHS > i_cut & PBS_value > p_cut)

## Write table with filters
write.table(x = final2, file = "filtered_pbs_vs_ihs.tsv", sep = "\t", 
            col.names = T, row.names = F, quote = F)

## Get a list with chromosomes for pbs
chromosomes <- significative_pbs %>%
  pull(CHROM) %>% 
  unique() %>% 
  sort()

## Get a list with chromosomes for ihs
chromosomes2 <- significative_ihs %>%
  pull(CHROM) %>% 
  unique() %>% 
  sort()

## Get a list with chromosomes for merging
chromosomes3 <- final2 %>%
  pull(CHROM) %>% 
  unique() %>% 
  sort()

## Making function for POS pbs
pos_pbs <- function(X) {
    significative_pbs %>% 
    filter(CHROM == X) %>% 
    pull(POS)
}

## App. to all chrs
map(chromosomes, pos_pbs) -> pbs_pos

##  Making function for PBS value
values_pbs <- function(X) {
    significative_pbs %>% 
    filter(CHROM == X) %>% 
    pull(PBS_value)
}

## App. to all chrs
map(chromosomes, values_pbs) -> pbs_values

## Making function for iHS pos
pos_ihs <- function(X) {
  significative_ihs %>% 
    filter(CHROM == X) %>% 
    pull(POS)
}

## App. to all chrs
map(chromosomes2, pos_ihs) -> ihs_pos

##  Making function for ihs value
values_ihs <- function(X) {
  significative_ihs %>% 
    filter(CHROM == X) %>% 
    pull(Std_iHS)
}

## App. to all chrs
map(chromosomes2, values_ihs) -> ihs_values

## Making function for iHSvsPBS pos
pos_ivp <- function(X) {
    final2 %>% 
      filter(CHROM == X) %>% 
      pull(POS)} 


## App. to all chrs
map(chromosomes3, pos_ivp) -> ivp_pos

##  Making function for iHSvsPBS value
values_ivp <- function(X) {
  final2 %>% 
    filter(CHROM == X) %>% 
    pull(value)
}

## App. to all chrs
map(chromosomes3, values_ivp) -> ivp_values

##  Making function for iHSvsPBS label
labels_ivp <- function(X) {
  final2 %>% 
    filter(CHROM == X) %>% 
    pull(label)
}

## App. to all chrs
map(chromosomes3, labels_ivp) -> ivp_labels

## Name the lists for reproducibility
names(pbs_pos) <- chromosomes
names(pbs_values) <- chromosomes
names(ihs_pos) <- chromosomes2
names(ihs_values) <- chromosomes2
names(ivp_pos) <- chromosomes3
names(ivp_values) <- chromosomes3
names(ivp_labels) <- chromosomes3

## Begin track with min chromosome
begginer <- min(chromosomes)

## First track PBS

tracks = BioCircosSNPTrack('SNPTrack', chromosomes = begginer, positions = pbs_pos[[paste(begginer)]],
                           values = pbs_values[[paste(begginer)]], colors = "#abdda4", size = 1.2,  
                           maxRadius = 0.97, minRadius = 0.82)

## Track for PBS
for (chr in chromosomes) {
  tracks = tracks + BioCircosSNPTrack('SNPTrack', chromosomes = chr, positions = pbs_pos[[paste(chr)]],
                                      values = pbs_values[[paste(chr)]], colors = "#abdda4", 
                                      maxRadius = 0.97, minRadius = 0.82, size = 1.2)
}

## Track for iHS
for (chr in chromosomes2) {
  tracks = tracks + BioCircosSNPTrack("SNPTrack2",chromosomes = chr, positions = ihs_pos[[paste(chr)]],
                                      values = ihs_values[[paste(chr)]], colors = "#ff764c", size = 1.5, 
                                      maxRadius = 0.79, minRadius = 0.64) 
}

## Track for merged
if (is_empty(chromosomes3) == F) {
for (chr in chromosomes3) {
  tracks = tracks + BioCircosSNPTrack("SNPTrack3",chromosomes = chr, positions = ivp_pos[[paste(chr)]],
                                      values = ivp_values[[paste(chr)]], colors = "#d53e4f", size = 1.5, 
                                      maxRadius = 0.61, minRadius = 0.46, labels = ivp_labels[[paste(chr)]]) 
  }
}

## Background tracks

tracks = tracks + BioCircosBackgroundTrack("testBGtrack1", minRadius = 0.82, maxRadius = 0.97,
                                           borderColors = "#F7FEF7", fillColors = "#F7FEF7", borderSize = 0.6)


tracks = tracks + BioCircosBackgroundTrack("testBGtrack1", minRadius = 0.64, maxRadius = 0.79,
                                           borderColors = "#FFEEBB", fillColor = "#FFEEBB", borderSize = 0.6)

tracks = tracks + BioCircosBackgroundTrack("testBGtrack1", minRadius = 0.46, maxRadius = 0.61,
                                           borderColors = "#FFEEEE", fillColor = "#FFEEEE", borderSize = 0.6)

##############################################################################

circleplot <- BioCircos(tracks, displayGenomeBorder = F)

##
htmltools::save_html(html = circleplot, file = "circleplot.html")
