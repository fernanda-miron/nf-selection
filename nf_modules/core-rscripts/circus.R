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
p_cut <- args[3]
i_cut <- args[4]

## Up data frame 
pbs <- vroom(pbs_file)
ihs <- vroom(ihs_file)

## Filter dataframe to only keep snps with higher PBS value
significative_pbs <- pbs %>% 
  filter(PBS_value > p_cut)

## Filter dataframe to only keep snps with higher iHS value
significative_ihs <- ihs %>% 
  filter(Std_iHS > i_cut)

## Filter to get SNPs from iHS and PBS
## First, lets get the top 1% of PBS
arrange_PBS <- pbs[with(pbs, order(-PBS_value)),]

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

## Add column with value and label
final2 <- final %>% 
  mutate(value = Std_iHS + PBS_value) %>% 
  mutate(label = paste("PBS =", PBS_value, "iHS =", Std_iHS))

## Filter selection
final3 <- final2 %>% 
  select(CHROM, POS, label)

## Write table with filters
write.table(x = final3, file = "filtered_pbs_vs_ihs.tsv", sep = "\t", 
            col.names = T, row.names = F, quote = F)

## Get a list with chromosomes for pbs
chromosomes <- significative_pbs %>%
  pull(CHROM)

## Get a list with chromosomes for ihs
chromosomes2 <- significative_ihs %>%
  pull(CHROM)

## Get a list with chromosomes for merging
chromosomes3 <- final2 %>%
  pull(CHROM)

## Get a list with pos for pbs
positions <- significative_pbs %>%
  pull(POS)

## Get a list with chromosomes for ihs
positions2 <- significative_ihs %>%
  pull(POS)

## Get a list with chromosomes for merging
positions3 <- final2 %>%
  pull(POS)

## Get a list with VALUE for pbs
value <- significative_pbs %>%
  pull(PBS_value)

## Get a list with chromosomes for ihs
value2 <- significative_ihs %>%
  pull(Std_iHS)

## Get a list with chromosomes for merging
value3 <- final2 %>%
  pull(value)

## Get a list with chromosomes for merging
labels <- final2 %>%
  pull(label)

## Track PBS
tracks = BioCircosSNPTrack('SNPTrack', chromosomes = chromosomes, positions = positions,
                           values = value, colors = "#abdda4", size = 1.2,  
                           maxRadius = 0.97, minRadius = 0.82, labels = "PBS")

## Track for iHS
tracks = tracks + BioCircosSNPTrack("SNPTrack2",chromosomes = chromosomes2, positions = positions2,
                                      values = value2, colors = "#ff764c", size = 1.2, 
                                      maxRadius = 0.79, minRadius = 0.64, labels = "iHS") 

## Track for merged
if (is_empty(chromosomes3) == F) {
  tracks = tracks + BioCircosSNPTrack("SNPTrack3",chromosomes = chromosomes3, positions = positions3,
                                      values = value3, colors = "#d53e4f", size = 3, 
                                      maxRadius = 0.61, minRadius = 0.46, labels = labels, shape = "rect") 
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
htmlwidgets::saveWidget(widget = circleplot, file = "circleplot.html")
