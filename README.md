# **nf_selection**


[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A521.10.5-green.svg)](https://www.nextflow.io/)
[![Ubuntu](https://img.shields.io/badge/ubuntu-%E2%89%A518.04-orange.svg)](https://ubuntu.com/download)
[![R](https://img.shields.io/badge/R-%E2%89%A54.1.2-blue.svg)](https://ubuntu.com/download)


### **Introduction**

**nf-selection** is a bioinformatics best-practice analysis pipeline for PBS (Population Branch Statistic)
and iHS (integrated Haplotype Score) computing to detect signals of natural selection. 

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner.

### **Pipeline Summary**

By default, the pipeline currently performs the following:

1. VCF Phasing (`SHAPEIT4`)
2. Ancestral annotation (`Jvarkit`)
3. iHS calculation (`rehh`)
4. iHS plotting (`ihs_rehh.R`)
5. Fst calculation (`vcftools`)
6. PBS calculation (`pbs_calculator.R`)
7. PBS and iHS cross-reference (`merging_pbs_ihs.R`)
8. Results annotation (`annovar`)

### **Requirements**


#### Compatible OS:



-   [Ubuntu 18.04 ](http://releases.ubuntu.com/18.04/)
-   [Ubuntu 20.04 ](http://releases.ubuntu.com/20.04/)



#### Software:



|                    Requirement                     |          Version           |  Required Commands \*  |
|:--------------------------------------------------:|:--------------------------:|:----------------------:|
|        [Nextflow](https://www.nextflow.io/)        |          21.10.5           |        nextflow        |
|          [R](https://www.r-project.org/)           |           4.1.2            |   PBS calculation, plotting, data wrangling   |
| [VCFtools](http://vcftools.sourceforge.net/)       |           0.1.15           | Fst calculation   |
| [bcftools](https://samtools.github.io/bcftools/bcftools.html)       |           1.10.2           | VCF ancestral annotation   |
| [samtools](http://www.htslib.org/)       |           1.12           | VCF ancestral annotation   |
|          [SHAPEIT4](https://odelaneau.github.io/shapeit4/)           |           4.2            | Haplotype phasing|


\* **Nextflow, R, VCFTools, and SHAPEIT4** must be accessible from your `$PATH` (*i.e.* you
should be able to use them in any path).

#### Databases:
On the following path "nf_modules\annovar", the user must download different databases as follows:

```
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
annotate_variation.pl -buildver hg38 -downdb cytoBand humandb/
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar exac03 humandb/ 
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp147 humandb/
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp30a humandb/
```

\* Due to license restrictions, users must download annotate_variation.pl at https://annovar.openbioinformatics.org/en/latest/user-guide/startup. User must also ensure the downloads are made on the correct genome reference. For more information check: https://annovar.openbioinformatics.org/en/latest/user-guide/startup/ (table_annovar.pl section). 

#### R packages needed:

|                    Requirement                     |          Use in workflow   |  
|:--------------------------------------------------:|:--------------------------:|
|        [stringr](https://cran.r-project.org/web/packages/stringr/index.html)        | Regular expression manipulation |
|          [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)          |   PBS calculation; plots development |
|          [tidyr](https://cran.r-project.org/web/packages/tidyr/index.html)          |   preprocessing and PBS calculation; data frame manipulation |
|          [cowplot](https://cran.r-project.org/web/packages/cowplot/index.html)          |   PBS calculation; plots development |
|        [rehh](https://cran.r-project.org/web/packages/rehh/index.html)        | iHS computation |
|        [vcfR](https://cran.r-project.org/web/packages/vcfR/vignettes/intro_to_vcfR.html)        | vcf manipulation |
|        [qqman](https://cran.r-project.org/web/packages/qqman/vignettes/qqman.html)        | Results visualization |
|        [scales](https://scales.r-lib.org/)        | Results visualization |
|        [vroom](https://www.tidyverse.org/blog/2019/05/vroom-1-0-0/)        | Dataframe reading |
|        [circlize](https://scales.r-lib.org/)        | Results visualization |

### **Installation**

Download nf_selection from Github repository:

    git clone https://github.com/fernanda-miron/nf-selection.git

------------------------------------------------------------------------

### **Test**
Before testing, user must unzip the following files: 

    gzip -d /test/data/ihs_files/chr2_1pops_lct.vcf.gz

    gzip -d /test/data/pbs_files/chr2_3pops_lct.vcf.gz
    
    gzip -d nf_modules/ancestral/ANCESTOR_for_chromosome_GRCh38_2_1_242193529_1.fa.gz



To test nf_selection execution using test data, run:

    bash runtest.sh

Your console should print the Nextflow log for the run, once every
process has been submitted, the following message will appear:

     ======
     Basic pipeline TEST SUCCESSFUL
     ======

results for test data should be in the following file:

    nf_selection/test/results

------------------------------------------------------------------------

### **Usage**

### **Main arguments**

You will need to create two design files with information about the files required 
for PBS and iHS respectively. You will need to use the following parameters to specify 
its location.

    --input_ihs [path to ihs design file]
    --input_pbs [path to pbs design file]

Both files have to be a comma-separated file as shown in the examples below. 

### **iHS design file**

The pipeline is designed to receive two different types of iHS design file 
depending on different situations:

If your VCF file is **NOT PHASED**, the design file should be as the example below:

    chromosome,path_vcf,path_genetic_map, path_reference_vcf, path_index_reference, path_ancestral, path_manifest 
    1,chr1.CEU.recode.vcf,genetic_map_chr1_combined_b37.txt, reference.vcf, reference.vcf.gz.tbi

|                    Column                     |          Description   |  
|:--------------------------------------------------:|:--------------------------:|
| chromosome | Chromosome number |
| path_vcf | Full path to VCF file per chromosomes. File MUST be unzipped |
| path_genetic_map | Full path to genetic_map.txt file. File MUST be unzipped |
| path_reference_vcf | Full path to reference VCF. File MUST be unzipped |
| path_index_reference | Full path to index file for reference VCF. File MUST be unzipped |
| path_ancestral | Full path to ancestral fasta. File MUST be unzipped. |
| path_manifest | Full path to manifest file for ancestral annotation. |


\* Ancestral files for GRCh38 can be found at nf-modules/ancestral

\* Manifest files for GRCh38 can be found at nf-modules/manifest 

----------------------------------------------------------------------------------------------------------

If your VCF file is **PHASED**, the design file should be as the example below:

    chromosome,path_vcf,path_genetic_map,path_ancestral,path_manifest
    2,test/data/ihs_files/ceu_filtered.recode.vcf,test/data/ihs_files/chr2.b38.predicted.map,nf_modules/ancestral_fasta/ANCESTOR_for_chromosome_GRCh38_2_1_242193529_1.fa,nf_modules/manifest_annotation/manifest2
    


\* The .map file must be the same as the one used for phasing

\* **As iHS is computed per chromosome, the user must provide the required info per chromosome**

----------------------------------------------------------------------------------------------------------

### **PBS design file**

The design file for PBS should be as the example below:

    path_vcf,path_pop1,path_pop2,path_popout
    all.chr21.recode.vcf,pop_1,pop_2,pop_out

|                    Column                     |          Description   |  
|:--------------------------------------------------:|:--------------------------:|
| path_vcf | Full path to VCF file. File MUST be unzipped |
| path_pop1 | Full path to file with the individuals IDs of the pop of interest |
| path_pop2 | Full path to file with the individuals IDs of the close related population|
| path_pop3 | Full path to file with the individuals IDs of the far related population|

\* As PBS is **NOT** computed per chromosome, the user must provide a VCF with merged chromosomes. 

\* The files with the individuals IDs MUST keep the names pop1, pop2 and pop3

------------------------------------------------------------------------

For information about options and parameters, run:

    nextflow run main.nf --help

------------------------------------------------------------------------

### **Other arguments**

Other arguments that the user may provide:

    --notphased [This argument may be used if the input VCF file is not phased]
    --notancestral [This argument may be used if the input VCF file is not ancestrally annotated]
    --genetic_map [This argument may be used if the user want to set its own genetic map]
    --imerged AND --pmerged [This arguments may be used TOGETHER to set the cutoff for iHS and PBS in the process of merging results]
    
-------------------------------------------------------------------------------
### **Running the pipeline**
The typical command for running the pipeline is as follows:

    nextflow run main.nf --input_ihs '[path to design file]' --input_pbs '[path to design file]'

### **Authors**

Fernanda Miron T

Israel Aguilar Ordonez