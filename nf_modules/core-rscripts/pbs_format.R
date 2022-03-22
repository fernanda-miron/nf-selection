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
#args[1] <- "./pbs.tsv"
#args[2] <- "./mart_export (2).txt"

## Place args into named object
file_pbs <- args[1]
file_biomart <- args[2]

## Upload files 
pbs.df <- vroom(file_pbs)
biomart.df <- vroom(file_biomart)

## Lets beggin with PBS file
## Change headers
names(pbs.df)[1] <- "seqname"
names(pbs.df)[2] <- "start"

## Add columns for GFF format
pbs.gff<- pbs.df %>% mutate(source = "PBS") %>%
  mutate(end = start) %>% 
  mutate(feature = "SNV") %>% 
  mutate(score = ".") %>% 
  mutate(strand = ".") %>% 
  mutate(frame = ".") %>% 
  mutate(attribute = str_c("PBS=", PBS_value))

## Modify order to match GFF order
pbs.gff<- pbs.gff %>% select(seqname, source, feature, start, end, score, 
                           strand, frame, attribute)

## Write GFF file
write.table(pbs.gff, 
            file = "pbs.gff", 
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
