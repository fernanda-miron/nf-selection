# This script can take Fst Data from VCFTools
# and calculate PBS values

# Charging library
library("dplyr")
library("ggplot2")
library("scales")
library("tidyr")
library("cowplot")

## Read args from command line
args <- commandArgs(trailingOnly = T)

# Uncomment for debbuging
# Comment for production mode only
#args[1] <- "."
#args[2] <- "pbs_snp.png"
#args[3] <- "pbs.tsv"

## Place args into named object
file_dir <- args[1]
jpg_file <- args[2]
tsv_file <- args[3]

## Using function for reading
fst_reader_snp<-function(filename,pops){
  fst<-read.table(file = filename, header = T, stringsAsFactors = F, sep = "\t") #read in fst output from vcftools
  colnames(fst)[3]<-paste(pops,".Fst",sep = "") #change fst column name to identify pops
  fst[,3][which(is.na(fst[,3]))]<-0 #change NA's to 0, the NA's are produced by vcftools when there is no variation at a site
  assign(pops,fst,envir = .GlobalEnv) #export dataframe with pop names as variable names
}

## Reading data
snp.name <- list.files(path = file_dir, pattern = "*.fst", full.names = T)

## Preparing data
pop1_snp <- as.list(strsplit(snp.name[1], "/") [[1]])
pop1s_saved <- as.character(pop1_snp[[2]])

pop2_snp <- as.list(strsplit(snp.name[2], "/") [[1]])
pop2s_saved <- as.character(pop2_snp[[2]])

pop3_snp <- as.list(strsplit(snp.name[3], "/") [[1]])
pop3s_saved <- as.character(pop3_snp[[2]])

## Reading data
psnp1 <- fst_reader_snp(filename = snp.name[1], pops = pop1s_saved)
psnp2 <- fst_reader_snp(filename = snp.name[2], pops = pop2s_saved)
psnp3 <- fst_reader_snp(filename = snp.name[3], pops = pop3s_saved)

# Make a list of your dataframes and join them all together
predata <- left_join(x = psnp1,
                     y = psnp2,
                     by = c("CHROM", "POS"))

all_fst <- left_join(x = predata,
                     y = psnp3,
                     by = c("CHROM", "POS"))

## PBS function
dopbs<-function(pop1_out,pop1_pop2,pop2_out) {
  Tpop1_out= -log(1-pop1_out)
  Tpop1_pop2= -log(1-pop1_pop2)
  Tpop2_out= -log(1-pop2_out)
  pbs= (Tpop1_out + Tpop1_pop2 - Tpop2_out)/2
  pbs
}

# Running for my data
my_pbs <-dopbs(pop1_out = all_fst[4],
              pop1_pop2 = all_fst[3],
              pop2_out = all_fst[5])

# Turn PBS into data frame and merge
pbsresults<-all_fst[,1:2]
pbsresults<-cbind(pbsresults,my_pbs)
pbsresults <- pbsresults %>% 
  rename(PBS_value = paste0(pop2s_saved,".Fst"))

#set negative PBS values to 0 for convenience
pbsresults$PBS_value[which(pbsresults$PBS_value<0)]<-0

# Arranging data to see values
arreglado <- pbsresults[order(-pbsresults$PBS_value),]

## Saving df
write.table(arreglado, file = tsv_file, col.names = T, row.names = F, sep = "\t")

## Making manhattan plot
## only 1-22, X o Y
valid_chroms <- c(1:22, "X", "Y")

#cleaning for only valid chrs 
clean_data.df <- arreglado %>%
  filter( CHROM %in% valid_chroms ) %>%
  ##Para poder ver chrom X y Y, los renumeramos a 23 y 24 respectivamente
  mutate( CHROM = gsub( pattern = "X", replacement = "23", x = CHROM ),
          CHROM = gsub( pattern = "Y", replacement = "24", x = CHROM )) %>%
  ## convertir en numericos los valores de chromosoma y posicion
  mutate( CHR_ID = as.numeric(CHROM),
          CHR_POS = as.numeric(POS))

##Preparamos la data para graficar en X un solo eje,
# a pesar de que sean diferentes cromosomas
single_axis.df <- clean_data.df %>% 
  # Compute chromosome size
  group_by( CHR_ID ) %>% 
  summarize(chr_len = max( CHR_POS) ) %>%
  # Calculate cumulative position of each chromosome
  mutate( chr_start_on_x = cumsum(chr_len) - chr_len) %>%
  select( -chr_len )

## Juntar clean data y el single axis para homogenizar la posicion
adjusted.df <- left_join( x = clean_data.df, y = single_axis.df,  by = "CHR_ID" ) %>%
  # ajustar posicion de acuerdo a un solo eje X
  arrange( CHR_ID , CHR_POS) %>%
  mutate( position_on_x = CHR_POS + chr_start_on_x,
          CHR_ID = as.factor(CHR_ID) )

## Vamos a crear un vector con las etiquetas para el eje X
xaxis.df = adjusted.df %>% group_by( CHR_ID ) %>%
  summarize( center = ( max( position_on_x ) + min( position_on_x) ) / 2 ) #Calculamos la posición central para cada cromosoma

##Limpiar dataframes que no se usaran
rm( clean_data.df )
rm( single_axis.df )

## Creating
man_uno.p <- ggplot( data = adjusted.df,
                     mapping = aes( x = position_on_x,
                                    y = PBS_value,     # transformamos los valores de Y a menos log10
                                    color = CHR_ID ) )   +
  geom_point( alpha = 0.1 )

## Creating basic plot
# Vamos a cambiar la escala de colores, a algo mas estandar
man_dos.p <- man_uno.p +
  scale_color_manual( values = rep(c("#F95738","#1B998B"), 12 ))

# vamos a poner los nombres de los cromosomas
man_tres.p <- man_dos.p +
  scale_x_continuous( label = xaxis.df$CHR_ID,      # Le ponemos un eje especial a X
                      breaks= xaxis.df$center,
                      expand = c(0.02, 0.5)) 

# Agregamos una linea roja con el corte de PBS
PBS_cutoff <- 0.8
man_cuatro.p <- man_tres.p +
  geom_hline( yintercept = PBS_cutoff,     
              color = "#2B2C28", lty = "dashed")        

# Ponemos titulo, etiqueta en X y modificamos tema
man_cinco.p <- man_cuatro.p +
  labs(title = "PBS by SNP´s",    # ponemos titulo dinamico con la funcion paste
       x = "CHR",
       y = "PBS value") +
  theme_light(base_size = 12) +                              
  theme( legend.position="none",            
         panel.grid = element_blank(),
         panel.spacing = ,)

# Separamos titulos de eje de plot
man_seis.p <- man_cinco.p +
  theme(axis.title.x = element_text(margin=margin(10,0,0,0), face = "bold", color = "grey20"),
        axis.title.y = element_text(margin=margin(0,10,0,0), face = "bold", color = "grey20"),
        plot.title=element_text(size=15,face="bold", color = "grey20"))

# Makin bar plot
p1 <- ggplot(data = arreglado, mapping = aes(x = PBS_value)) +
  geom_histogram(color="#EE964B", fill="#EE964B", alpha=0.2) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(expand = c(0.02, 0)) +
  theme_bw() +
  theme(axis.title.x = element_text(margin=margin(10,0,0,0), face = "bold", color = "grey20"),
        axis.title.y = element_text(margin=margin(0,10,0,0), face = "bold", color = "grey20"),
        plot.title=element_text(size=15,face="bold", color = "grey20")) +
  labs(title = paste("PBS by SNP"),
       y = "Number of SNPs",
       x = "PBS value") 
p1

# Spiderplot
## Making spiderplot graph
## Reading data
frq.names <- list.files(path = file_dir, pattern = "*.frq", full.names = T)

## Reading first data frame
AF_1 <- read.table(file = frq.names[1], sep = "\t",
                   header = T, stringsAsFactors = F, 
                   fill = T, row.names = NULL)

## Changing df format
colnames(AF_1) <- c("CHROM", "POS", "N_ALLELES", "N_CHR", "AF1_1", "AF2_1")
AF_1 <- transform(AF_1, CHROM = as.numeric(CHROM))
AF_1 <- transform(AF_1, POS = as.numeric(POS))

## Reading second dataframe
AF_2 <- read.table(file = frq.names[2], sep = "\t",
                   header = T, stringsAsFactors = F, 
                   fill = T, row.names = NULL)

## Changing df format
colnames(AF_2) <- c("CHROM", "POS", "N_ALLELES", "N_CHR", "AF1_2", "AF2_2")
AF_2 <- transform(AF_2, CHROM = as.numeric(CHROM))
AF_2 <- transform(AF_2, POS = as.numeric(POS))

## Reading third dataframe
AF_3 <- read.table(file = frq.names[3], sep = "\t",
                   header = T, stringsAsFactors = F, 
                   fill = T, row.names = NULL)

## Changing df format
colnames(AF_3) <- c("CHROM", "POS", "N_ALLELES", "N_CHR", "AF1_3", "AF2_3")
AF_3 <- transform(AF_3, CHROM = as.numeric(CHROM))
AF_3 <- transform(AF_3, POS = as.numeric(POS))

## Merge data
merged_data.df <- arreglado %>% left_join(AF_1,
                                          by = c("CHROM"="CHROM", "POS"="POS")) %>%
  left_join(AF_2, by = c("CHROM"="CHROM","POS"="POS","N_ALLELES" = "N_ALLELES")) %>%
  left_join(AF_3, by = c("CHROM"="CHROM","POS"="POS","N_ALLELES" = "N_ALLELES")) %>%
  select(2:3,6:7,9:10,12:13)

## Changing formats
fixed_data.df <- merged_data.df %>%
  mutate(AF1_1 = gsub(x = AF1_1,
                      pattern = ".*:",
                      replacement = "")) %>%
  mutate(AF2_1 = gsub(x = AF2_1,
                      pattern = ".*:",
                      replacement = "")) %>%
  mutate(AF1_2 = gsub(x = AF1_2,
                      pattern = ".*:",
                      replacement = "")) %>%
  mutate(AF2_2 = gsub(x = AF2_2,
                      pattern = ".*:",
                      replacement = "")) %>%
  mutate(AF1_3 = gsub(x = AF1_3,
                      pattern = ".*:",
                      replacement = "")) %>%
  mutate(AF2_3 = gsub(x = AF2_3,
                      pattern = ".*:",
                      replacement = ""))

## Add SNP´s name
nombres <- sprintf("SNP%s",seq(1:nrow(fixed_data.df)))
fixed_data.df$SNP <- nombres

## Pivot dataframe
fixed_data.df <- pivot_longer(fixed_data.df, cols = c("AF1_1", "AF1_2", "AF1_3",
                                                      "AF2_1", "AF2_2", "AF2_3" ),
                              names_to = "AF", values_to = "valor")
## FIltering
row_name <- "SNP1"
fst_values <- c("AF2_1", "AF2_2", "AF2_3")
fixed_data.df <- transform(fixed_data.df, valor = as.numeric(valor))

## Ploting first spider plot
spider_uno.p <- fixed_data.df %>%
  filter(SNP == row_name) %>%
  filter(AF == fst_values) %>%
  ggplot( mapping = aes(x = AF, y = valor) ) +
  geom_point( size = 3, color = "#49697F" )
spider_uno.p

## Using geom segment geometry
spider_dos.p <- spider_uno.p +
  geom_segment(
    aes( x = AF, xend = AF,
         y = 0.0, yend = valor), color = "#49697F", size = 1
  )
spider_dos.p

## Coord polar
spider_tres.p <- spider_dos.p + 
  coord_polar()
spider_tres.p

## Adding value
spider_cuatro.p <- spider_tres.p +
  geom_text( aes(label = valor))
spider_cuatro.p

## Improving geometry
spider_cuatro.p <- spider_tres.p +
  geom_text( aes(label = valor),position = position_nudge(y = 0.3))
spider_cuatro.p

# Cleaning plot
spider_cinco.p <- spider_cuatro.p +
  theme_light() +                         
  theme(panel.border = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank(),       
        axis.text.y = element_blank(),       
        axis.title = element_blank()         
  )
spider_cinco.p

## Adding dinamic title
spider_seis.p <- spider_cinco.p +
  labs(title = "SNP with higher PBS value") +
  xlab("Number of SNPs by 10KB")+
  ylab("PBS") +
  theme(plot.title=element_text(size=15,face="bold", color = "grey20"))
spider_seis.p

## Merging
grid1 <- plot_grid(p1, spider_seis.p, align = "h", labels = c('B', 'C'))
grid2 <- plot_grid(man_seis.p, grid1, align = "h", labels = 'A', nrow = 2)

# Guardar plot
ggsave(filename = jpg_file, 
       plot = grid2, 
       device = "png",
       width = 15, height = 11, units = "in",
       bg = "white",
       dpi = 300)
