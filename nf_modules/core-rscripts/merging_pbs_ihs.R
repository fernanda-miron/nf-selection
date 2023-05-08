## This script takes ihs and pbs data and makes a circus plot
library("dplyr")
library("vroom")
library("circlize")

## Read args from command line
args <- commandArgs(trailingOnly = T)

# Uncomment for debbuging
# Comment for production mode only
#args[1] <- "one_percent_pbs.tsv"
#args[2] <- "final_ihs_onepercent.tsv"

## Place args into named object
pbs_file <- args[1]
ihs_file <- args[2]

## Up data frame 
pbs <- vroom(pbs_file)
ihs <- vroom(ihs_file)

## Rename for merging
ihs <- ihs %>% 
  rename("CHROM" = CHR, "POS" = POSITION)

## Merge top PBS and top iHS
merged <- inner_join(pbs, ihs, by = c("POS", "CHROM"))

## Write merged table
write.table(x = merged, file = "pbs_vs_ihs.tsv", sep = "\t", 
            col.names = T, row.names = F, quote = F)

## Making bed
bed1 <- merged %>% 
  select(CHROM, POS)

## Adding same column
bed2 <- bed1 %>% 
  mutate(POS2 = POS)

## Save file
write.table(x = bed2, file = "final_bed", sep = "\t", 
            col.names = F, row.names = F, quote = F)

## Treating 
## First rename
final_pbs <- pbs %>% 
  rename("chr" = CHROM, "start" = POS, "value1" = PBS_value)

## Add end
final_pbs <- final_pbs %>% 
  mutate("end" = start)

## Select
final_pbs <- final_pbs %>% 
  select(chr, start, end, value1)

## Add chr to beginning
final_pbs$chr <- sub("^", "chr", final_pbs$chr)

## Treating 
## First rename
final_ihs <- ihs %>% 
  rename("chr" = CHROM, "start" = POS, "value1" = IHS)

## Add end
final_ihs <- final_ihs %>% 
  mutate("end" = start)

## Select
final_ihs <- final_ihs %>% 
  select(chr, start, end, value1)

## Add chr to beginning
final_ihs$chr <- sub("^", "chr", final_ihs$chr)

## 
png(filename="circus.png", width = 200, height = 200,
    units='mm', res = 1000, bg = "white")

###
circos.par("start.degree" = 90)
circos.initializeWithIdeogram(species = "hg38", chromosome.index = paste0("chr", c(1:22)))

## 
circos.genomicTrack(final_pbs, panel.fun = function(region, value, ...) {
  circos.genomicPoints(region, value, pch = 16, cex = 0.5, col = "#91C4F2", ...)
})


circos.genomicTrack(final_ihs, panel.fun = function(region, value, ...) {
  circos.genomicPoints(region, value, pch = 16, cex = 0.5, col = "#ACB2B5", ...)
})

##
dev.off()

## I should do this to clean plots
circos.clear()
