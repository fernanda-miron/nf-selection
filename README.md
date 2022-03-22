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

1. VCF Phasing (`SHAPEIT2`)
2. Ancestral annotation (`annotate.py`)
3. Map file generation (`make_map.py`)
4. iHS calculation (`hapbin`)
5. iHS plotting (`ihs_treatment.R`)
6. iHS annotation (`bedtools`)
7. Fst calculation (`vcftools`)
8. PBS calculation (`pbs_calculator.R`)
9. PBS annotation (`bedtools`)
10. PBS and iHS comparison (`R script - To do`)

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
|          [SHAPEIT2](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html)           |           2.17            | Haplotype phasing|
| [annotate.py](https://github.com/MerrimanLab/selectionTools/tree/master/selection_pipeline)       |           1.1           | Ancestral genome annotation   |
|        [make_map.py](https://github.com/evotools/hapbin/blob/master/tools/make_map.py)        |          NA           | Building genetic map        |
|          [hapbin](https://github.com/evotools/hapbin)           |           NA            | iHS computing  |



\* These commands must be accessible from your `$PATH` (*i.e.* you
should be able to invoke them from your command line).



#### R packages needed:

|                    Requirement                     |          Use in workflow   |  
|:--------------------------------------------------:|:--------------------------:|
|        [stringr](https://cran.r-project.org/web/packages/stringr/index.html)        | Regular expression manipulation |
|          [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)          |   PBS calculation; data frame manipulation |
|          [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)          |   PBS calculation; plots development |
|          [ggrepel](https://cran.r-project.org/web/packages/ggrepel/index.html)          |   PBS calculation; plots development |
|          [tidyr](https://cran.r-project.org/web/packages/tidyr/index.html)          |   preprocessing and PBS calculation; data frame manipulation |
|          [cowplot](https://cran.r-project.org/web/packages/cowplot/index.html)          |   PBS calculation; plots development |
|          [vroom](https://cran.r-project.org/web/packages/vroom/index.html)          |   PBS calculation; import of data frames |




\* You need to install the dplyr, ggplot2, ggrepel, tidyr and cowplot libraries in R/RStudio


### **Installation**

Download nf_selection from Github repository:

    git clone git@github.com:fernanda-miron/nf_selection

------------------------------------------------------------------------

### **Test**

To test nf_selection execution using test data, run:

    ./runtest.sh

Your console should print the Nextflow log for the run, once every
process has been submitted, the following message will appear:

     ======
     Basic pipeline TEST SUCCESSFUL
     ======

results for test data should be in the following file:

    nf_selection/test/results

------------------------------------------------------------------------

### **Usage**

To run with your own data go to the pipeline directory and execute:

    nextflow run ${pipeline_name}.nf --input <path to design csv file> [--output_dir path to results ]

For information about options and parameters, run:

    nextflow run wrangling.nf --help

------------------------------------------------------------------------

### **Authors**

Fernanda Miron T
Israel Aguilar Ordonez

