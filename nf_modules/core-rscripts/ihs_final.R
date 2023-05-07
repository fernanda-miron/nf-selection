## This script prepares iHS data to do PBS comparison

## Charging library
library("vroom")
library("dplyr")
library("stringr")
library("ggplot2")
library("scales")
library("cowplot")
library("qqman")

## Read args from command line
args <- commandArgs(trailingOnly = T)

# Uncomment for debbuging
# Comment for production mode only
#args[1] <- "."
#args[2] <- "./final_ihs.tsv"
#args[3] <- "./final_one_percent.tsv"
#args[4] <- "./ihs_histogram.png"


## Place args into named object
file_dir <- args[1]
tsv_file <- args[2]
tsv_file_onepercent <- args[3]
histogram_file <- args[4]

## Get all the ihs_file files in path
temp = list.files(pattern="*.tsv")

## Read each file as dataframe
for (i in 1:length(temp)) assign(temp[i], vroom(temp[i], col_names = T))

## Change header of last column of each df
## First, get a list with all de dfs in environment
dfs <- Filter(function(x) is(x, "data.frame"), mget(ls()))

## Get modified dataframes (again, pff)
dfs_edited <- Filter(function(x) is(x, "data.frame"), mget(ls()))

## Finally, bind dataframes
all_ihs <- do.call(rbind, dfs_edited)

## Getting absolute value of iHS
all_ihs$IHS <- abs(all_ihs$IHS)

## Get one percent
ihs_onepercent_sites <-all_ihs[which(all_ihs$pvals<=quantile(all_ihs$pvals,0.01)),]

## Save table with whole data
write.table(x = all_ihs, file = tsv_file, 
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

## Save table with one percent
write.table(x = ihs_onepercent_sites, file = tsv_file_onepercent, 
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

## 
## plot p -values
png("manhattan.png", width = 365, height = 225, units='mm', res = 300)

p1 <- manhattan(x = all_ihs, chr = "CHR", bp = "POSITION", p = "pvals", 
                snp = "POSITION", ylim = c(0, 10), main = "iHS values distribution", 
                col = c("#919098", "#91C4F2"), suggestiveline = F, genomewideline = F)
p1

dev.off()


# Making bar plot
p1 <- ggplot(data = all_ihs, mapping = aes(x = IHS)) +
  geom_histogram( color="#919098", fill="#919098", alpha=0.2) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(expand = c(0.02, 0)) +
  theme_light(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(margin=margin(10,0,0,0), face = "bold", color = "grey20"),
        axis.title.y = element_text(margin=margin(0,10,0,0), face = "bold", color = "grey20"),
        plot.title=element_text(size=15,face="bold", color = "grey20")) +
  labs(title = paste("iHS values"),
       y = "Number of SNPs",
       x = "|iHS|")
p1

# Finally, histogram plot
ggsave(filename = histogram_file, 
       plot = p1, 
       device = "png",
       width = 8, height = 6, units = "in",
       bg = "white",
       dpi = 300)
