###############################
  #CHANGE DATA TO GFF FORMAT#
###############################

#change PBS.df format to GFF
#Fields must be tab-separated. Also, all but the final field in each feature line must contain a value; "empty" columns should be denoted with a '.'
#seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. 
#source - name of the program that generated this feature, or the data source (database or project name)
#feature - feature type name, e.g. Gene, Variation, Similarity
#start - Start position* of the feature, with sequence numbering starting at 1.
#end - End position* of the feature, with sequence numbering starting at 1.
#score - A floating point value.
#strand - defined as + (forward) or - (reverse).
#frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
#attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.

## Charging libraries
library("dplyr")
library("vroom")
library("stringr")

## Read args from command line
args <- commandArgs(trailingOnly = T)

# Uncomment for debbuging
# Comment for production mode only
#args[1] <- "./final_ihs.tsv"
#args[2] <- "./mart_export.txt"

## Place args into named object
file_ihs <- args[1]
file_biomart <- args[2]

## Upload files 
ihs.df <- vroom(file_ihs)
biomart.df <- vroom(file_biomart)

## Select columns
sampled.df <- ihs.df %>% 
  select(chr, POS, Std_iHS)

## Lets beggin with PBS file
## Change headers
names(sampled.df)[1] <- "seqname"
names(sampled.df)[2] <- "start"

## Add columns for GFF format
ihs.gff<- sampled.df %>% mutate(source = "iHS") %>%
  mutate(end = start) %>% 
  mutate(feature = "SNV") %>% 
  mutate(score = ".") %>% 
  mutate(strand = ".") %>% 
  mutate(frame = ".") %>% 
  mutate(attribute = str_c("iHS=", Std_iHS))

## Modify order to match GFF order
ihs.gff<- ihs.gff %>% select(seqname, source, feature, start, end, score, 
                           strand, frame, attribute)

## Write GFF file
write.table(ihs.gff, 
            file = "ihs.gff", 
            quote = F,
            sep = "\t", 
            row.names = F,
            col.names = F)

## Give biomart.df GFF format
## First, rename columns
names(biomart.df)[1] <- "seqname"
names(biomart.df)[2] <- "gene_ID"
names(biomart.df)[3] <- "start"
names(biomart.df)[4] <- "end"
names(biomart.df)[5] <- "name"

## Add columns for GFF format
biomart.gff <- biomart.df %>%  mutate(source = "biomart") %>% 
  mutate(feature = "gene") %>% 
  mutate(score = ".") %>% 
  mutate(strand = ".") %>% 
  mutate(frame = ".") %>% 
  mutate(attribute = str_c("gene_id \"",gene_ID,"\"; gene_name \"", name, "\""))

## Change order to match the one of GFF
biomart.gff <- biomart.gff %>% select(seqname, source, feature, start, end, score, 
                                   strand, frame, attribute)

## Write GFF from mart
write.table(biomart.gff, 
            file = "biomart.gff", 
            quote = F,
            sep = "\t", 
            row.names = F,
            col.names = F)
