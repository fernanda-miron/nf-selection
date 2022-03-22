## This script prepares iHS data to do PBS comparison

## Charging library
library("vroom")
library("dplyr")
library("stringr")
library("ggplot2")
library("scales")
library("cowplot")

## Read args from command line
args <- commandArgs(trailingOnly = T)

# Uncomment for debbuging
# Comment for production mode only
#args[1] <- "."
#args[2] <- "./final_ihs.tsv"
#args[3] <- "2"
#args[4] <- "./ihs_manhattan.png"
#args[5] <- "./ihs_histogram.png"

## Place args into named object
file_dir <- args[1]
tsv_file <- args[2]
ihs_cut <- args[3]
png_file <- args[4]
histogram_file <- args[5]

## Get all the ihs_file files in path
temp = list.files(pattern="*.ihs_file")

## Read each file as dataframe
for (i in 1:length(temp)) assign(temp[i], vroom(temp[i], col_names = T))

## Change header of last column of each df
## First, get a list with all de dfs in environment
dfs <- Filter(function(x) is(x, "data.frame"), mget(ls()))

## Second, define col names
colnames <- c("Index","ID","Freq","iHH_0","iHH_1","iHS","Std_iHS","chr")

## Third, operate list to change col names. Ready!
list2env(lapply(dfs, setNames, colnames), .GlobalEnv)

## Get modified dataframes (again, pff)
dfs_edited <- Filter(function(x) is(x, "data.frame"), mget(ls()))

## Finally, bind dataframes
all_ihs <- do.call(rbind, dfs_edited)

## Get position for proper manage
## getting position from ID
changed_vcf <- str_split_fixed(all_ihs$ID, ":", 4)

## changing from matrix to dataframe format 
changed_vcf2 <- as.data.frame(changed_vcf) 

## Getting the column of interest
colum_of_interest <- changed_vcf2 %>% 
  select(V2) %>% 
  rename(POS = V2)

## Merging column of interest to original data 
final_ihs_vcf <- all_ihs %>% 
  mutate(POS = colum_of_interest$POS)

## Getting absolute value of iHS
final_ihs_vcf$Std_iHS <- abs(final_ihs_vcf$Std_iHS)
final_ihs_vcf$iHS <- abs(final_ihs_vcf$iHS)

## Lets save the table
write.table(x = final_ihs_vcf, file = tsv_file, 
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

## Lets beggin plotting
## First, a manhattan plot
## only 1-22, X o Y
valid_chroms <- c(1:22, "X", "Y")

#cleaning for only valid chrs 
clean_data.df <- final_ihs_vcf %>%
  filter( chr %in% valid_chroms ) %>%
  ## Rename X and Y chr as 23 and 24 respect.
  mutate( chr = gsub( pattern = "X", replacement = "23", x = chr ),
          chr = gsub( pattern = "Y", replacement = "24", x = chr )) %>%
  ## change pos and chrid as numeric
  mutate( CHR_ID = as.numeric(chr),
          CHR_POS = as.numeric(POS))

## Prepare data for single X axis plotting
# despite different chromosomes
single_axis.df <- clean_data.df %>% 
  # Compute chromosome size
  group_by( CHR_ID ) %>% 
  summarize(chr_len = max( CHR_POS) ) %>%
  # Calculate cumulative position of each chromosome
  mutate( chr_start_on_x = cumsum(chr_len) - chr_len) %>%
  select( -chr_len )

## Merge clean data and single axis
adjusted.df <- left_join( x = clean_data.df, y = single_axis.df,  by = "CHR_ID" ) %>%
  arrange( CHR_ID , CHR_POS) %>%
  mutate( position_on_x = CHR_POS + chr_start_on_x,
          CHR_ID = as.factor(CHR_ID) )

## Creating a vector with x axis labbels
xaxis.df = adjusted.df %>% group_by( CHR_ID ) %>%
  summarize( center = ( max( position_on_x ) + min( position_on_x) ) / 2 ) 

## Remove useless dfs
rm( clean_data.df )
rm( single_axis.df )

## Creating
man_uno.p <- ggplot( data = adjusted.df,
                     mapping = aes( x = position_on_x,
                                    y = Std_iHS,     
                                    color = CHR_ID ) )   +
  geom_point( alpha = 0.1 )

## Creating basic plot
# Changing color scale
man_dos.p <- man_uno.p +
  scale_color_manual( values = rep(c("#F95738","#1B998B"), 12 ))

# Adding chr names
man_tres.p <- man_dos.p +
  scale_x_continuous( label = xaxis.df$CHR_ID,      # Le ponemos un eje especial a X
                      breaks= xaxis.df$center,
                      expand = c(0.02, 0.5)) 

# Adding ihs cutoff
ihs_cutoff <- as.integer(ihs_cut)
man_cuatro.p <- man_tres.p +
  geom_hline( yintercept = ihs_cutoff,     
              color = "#2B2C28", lty = "dashed")  

# Adding title, X labbels and theme
man_cinco.p <- man_cuatro.p +
  labs(title = "iHS values",    
       x = "CHR",
       y = "|iHS|") +
  theme_light(base_size = 12) +                              
  theme( legend.position="none",            
         panel.grid = element_blank(),
         panel.spacing = ,)

# Making plot more beautiful
man_seis.p <- man_cinco.p +
  theme(axis.title.x = element_text(margin=margin(10,0,0,0), face = "bold", color = "grey20"),
        axis.title.y = element_text(margin=margin(0,10,0,0), face = "bold", color = "grey20"),
        plot.title=element_text(size=15,face="bold", color = "grey20"))

# Finally, Manhattan plot
ggsave(filename = png_file, 
       plot = man_seis.p, 
       device = "png",
       width = 15, height = 8, units = "in",
       bg = "white",
       dpi = 300)

# Making bar plot
p1 <- ggplot(data = final_ihs_vcf, mapping = aes(x = Std_iHS)) +
  geom_histogram( color="#EE964B", fill="#EE964B", alpha=0.2) +
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

# Finally, histogram plot
ggsave(filename = histogram_file, 
       plot = p1, 
       device = "png",
       width = 8, height = 6, units = "in",
       bg = "white",
       dpi = 300)
