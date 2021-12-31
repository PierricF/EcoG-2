CC2
================

# Installation des packages

``` bash
sudo apt-get update -y 
sudo apt-get install -y libbz2-dev
sudo apt-get install -y liblzma-dev
```

    ## sudo: unable to resolve host 012e23cc2553: Name or service not known
    ## Get:1 http://security.ubuntu.com/ubuntu focal-security InRelease [114 kB]
    ## Hit:2 http://archive.ubuntu.com/ubuntu focal InRelease
    ## Get:3 http://archive.ubuntu.com/ubuntu focal-updates InRelease [114 kB]
    ## Get:4 http://archive.ubuntu.com/ubuntu focal-backports InRelease [108 kB]
    ## Fetched 336 kB in 1s (394 kB/s)
    ## Reading package lists...
    ## sudo: unable to resolve host 012e23cc2553: Name or service not known
    ## Reading package lists...
    ## Building dependency tree...
    ## Reading state information...
    ## libbz2-dev is already the newest version (1.0.8-2).
    ## 0 upgraded, 0 newly installed, 0 to remove and 17 not upgraded.
    ## sudo: unable to resolve host 012e23cc2553: Name or service not known
    ## Reading package lists...
    ## Building dependency tree...
    ## Reading state information...
    ## liblzma-dev is already the newest version (5.2.4-1ubuntu1).
    ## 0 upgraded, 0 newly installed, 0 to remove and 17 not upgraded.

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("BiocStyle")
```

    ## 'getOption("repos")' replaces Bioconductor standard repositories, see
    ## '?repositories' for details
    ## 
    ## replacement repositories:
    ##     CRAN: https://packagemanager.rstudio.com/all/__linux__/focal/latest

    ## Bioconductor version 3.14 (BiocManager 1.30.16), R 4.1.2 (2021-11-01)

    ## Warning: package(s) not installed when version(s) same as current; use `force = TRUE` to
    ##   re-install: 'BiocStyle'

    ## Installation paths not writeable, unable to update packages
    ##   path: /usr/local/lib/R/library
    ##   packages:
    ##     Matrix

``` r
BiocManager::install("Rhtslib")
```

    ## 'getOption("repos")' replaces Bioconductor standard repositories, see
    ## '?repositories' for details
    ## 
    ## replacement repositories:
    ##     CRAN: https://packagemanager.rstudio.com/all/__linux__/focal/latest
    ## 
    ## Bioconductor version 3.14 (BiocManager 1.30.16), R 4.1.2 (2021-11-01)

    ## Warning: package(s) not installed when version(s) same as current; use `force = TRUE` to
    ##   re-install: 'Rhtslib'

    ## Installation paths not writeable, unable to update packages
    ##   path: /usr/local/lib/R/library
    ##   packages:
    ##     Matrix

``` r
library("knitr")
library("BiocStyle")
.cran_packages <- c("ggplot2", "gridExtra", "devtools")
install.packages(.cran_packages) 
```

    ## Installing packages into '/usr/local/lib/R/site-library'
    ## (as 'lib' is unspecified)

``` r
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn")
BiocManager::install(.bioc_packages)
```

    ## 'getOption("repos")' replaces Bioconductor standard repositories, see
    ## '?repositories' for details
    ## 
    ## replacement repositories:
    ##     CRAN: https://packagemanager.rstudio.com/all/__linux__/focal/latest

    ## Bioconductor version 3.14 (BiocManager 1.30.16), R 4.1.2 (2021-11-01)

    ## Warning: package(s) not installed when version(s) same as current; use `force = TRUE` to
    ##   re-install: 'dada2' 'phyloseq' 'DECIPHER' 'phangorn'

    ## Installation paths not writeable, unable to update packages
    ##   path: /usr/local/lib/R/library
    ##   packages:
    ##     Matrix

``` r
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)
```

    ##   ggplot2 gridExtra  devtools     dada2  phyloseq  DECIPHER  phangorn 
    ##      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE

``` r
install.packages("vegan")
```

    ## Installing package into '/usr/local/lib/R/site-library'
    ## (as 'lib' is unspecified)

``` r
install.packages("dplyr")
```

    ## Installing package into '/usr/local/lib/R/site-library'
    ## (as 'lib' is unspecified)

``` r
install.packages("venn")
```

    ## Installing package into '/usr/local/lib/R/site-library'
    ## (as 'lib' is unspecified)

# Chargement des packages

``` r
library("knitr")
library("BiocStyle")
.cran_packages <- c("ggplot2", "gridExtra")
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn")
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)
```

    ##   ggplot2 gridExtra     dada2  phyloseq  DECIPHER  phangorn 
    ##      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE

``` r
library("vegan")
library("dplyr")
library("venn")
```

``` r
set.seed(100)
```

# Préparation des données

On télécharge les données dans /EcoG-2/CC2_data.  
On indique le chemin vers nos données.

``` r
path <- "/home/rstudio/EcoG-2/CC2_data"
list.files(path)
```

    ##  [1] "filtered"               "SRR14295221"            "SRR14295221_1.fastq.gz"
    ##  [4] "SRR14295221_2.fastq.gz" "SRR14295222"            "SRR14295222_1.fastq.gz"
    ##  [7] "SRR14295222_2.fastq.gz" "SRR14295223"            "SRR14295223_1.fastq.gz"
    ## [10] "SRR14295223_2.fastq.gz" "SRR14295224"            "SRR14295224_1.fastq.gz"
    ## [13] "SRR14295224_2.fastq.gz" "SRR14295225"            "SRR14295225_1.fastq.gz"
    ## [16] "SRR14295225_2.fastq.gz" "SRR14295226"            "SRR14295226_1.fastq.gz"
    ## [19] "SRR14295226_2.fastq.gz" "SRR14295227"            "SRR14295227_1.fastq.gz"
    ## [22] "SRR14295227_2.fastq.gz" "SRR14295228"            "SRR14295228_1.fastq.gz"
    ## [25] "SRR14295228_2.fastq.gz" "SRR14295229"            "SRR14295229_1.fastq.gz"
    ## [28] "SRR14295229_2.fastq.gz" "SRR14295230"            "SRR14295230_1.fastq.gz"
    ## [31] "SRR14295230_2.fastq.gz" "SRR14295231"            "SRR14295231_1.fastq.gz"
    ## [34] "SRR14295231_2.fastq.gz" "SRR14295232"            "SRR14295232_1.fastq.gz"
    ## [37] "SRR14295232_2.fastq.gz" "SRR14295233"            "SRR14295233_1.fastq.gz"
    ## [40] "SRR14295233_2.fastq.gz" "SRR14295234"            "SRR14295234_1.fastq.gz"
    ## [43] "SRR14295234_2.fastq.gz" "SRR14295235"            "SRR14295235_1.fastq.gz"
    ## [46] "SRR14295235_2.fastq.gz"

Répartir les reads, extraire les noms,

``` r
fnFs <- sort(list.files(path, pattern = "_1.fastq.gz"))
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz"))

sampleNames <- sapply(strsplit(fnFs,"_"),'[',1)    
#boues <- c(sampleNames[6:7],sampleNames[9],sampleNames[13:14])
#déjections <- c(sampleNames[8],sampleNames[10:12],sampleNames[15])
#compost <- sampleNames[1:5]
#sampleNames <- c(boues,déjections,compost)

fnFs <- file.path(path,fnFs)  
fnRs <- file.path(path,fnRs)
```

Inspection des profils qualité

``` r
plotQualityProfile(fnFs[1:2])
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

![](CC2_files/figure-gfm/unnamed-chunk-8-1.png)<!-- --> On a une bonne
qualité pour les reads forward. On choisit d’ôter quelques derniers
nucléotides en coupant à la 140e position en plus d’enlever les 10
premières bases, qui sont plus susceptibles de contenir des erreurs.

``` r
plotQualityProfile(fnRs[1:2])
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

![](CC2_files/figure-gfm/unnamed-chunk-9-1.png)<!-- --> On a une bonne
qualité pour les reads reverse. On choisit d’ôter quelques derniers
nucléotides en coupant à la 130e position en plus d’enlever les 10
premières bases, qui sont plus susceptibles de contenir des erreurs.

Définir les noms de fichiers pour les fastq.gz filtrés

``` r
filt_path <- file.path(path,"filtered") 
if(!file_test("-d",filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path,paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path,paste0(sampleNames, "_R_filt.fastq.gz"))
```

Filtrer les reads

``` r
out <- filterAndTrim(fnFs,filtFs,fnRs,filtRs,truncLen=c(140,130),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)
```

    ##                        reads.in reads.out
    ## SRR14295221_1.fastq.gz     9819      9498
    ## SRR14295222_1.fastq.gz     5813      5674
    ## SRR14295223_1.fastq.gz    11669     11310
    ## SRR14295224_1.fastq.gz     7251      7029
    ## SRR14295225_1.fastq.gz     6606      6435
    ## SRR14295226_1.fastq.gz     7496      7301

# Définir les ASVs

Dérépliquer

``` r
derepFs <- derepFastq(filtFs, verbose=TRUE)  
```

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/EcoG-2/CC2_data/filtered/SRR14295221_F_filt.fastq.gz

    ## Encountered 3524 unique sequences from 9498 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/EcoG-2/CC2_data/filtered/SRR14295222_F_filt.fastq.gz

    ## Encountered 1784 unique sequences from 5674 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/EcoG-2/CC2_data/filtered/SRR14295223_F_filt.fastq.gz

    ## Encountered 3893 unique sequences from 11310 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/EcoG-2/CC2_data/filtered/SRR14295224_F_filt.fastq.gz

    ## Encountered 2261 unique sequences from 7029 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/EcoG-2/CC2_data/filtered/SRR14295225_F_filt.fastq.gz

    ## Encountered 2470 unique sequences from 6435 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/EcoG-2/CC2_data/filtered/SRR14295226_F_filt.fastq.gz

    ## Encountered 2609 unique sequences from 7301 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/EcoG-2/CC2_data/filtered/SRR14295227_F_filt.fastq.gz

    ## Encountered 2571 unique sequences from 7010 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/EcoG-2/CC2_data/filtered/SRR14295228_F_filt.fastq.gz

    ## Encountered 2755 unique sequences from 9493 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/EcoG-2/CC2_data/filtered/SRR14295229_F_filt.fastq.gz

    ## Encountered 2612 unique sequences from 7002 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/EcoG-2/CC2_data/filtered/SRR14295230_F_filt.fastq.gz

    ## Encountered 2758 unique sequences from 8647 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/EcoG-2/CC2_data/filtered/SRR14295231_F_filt.fastq.gz

    ## Encountered 2505 unique sequences from 7932 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/EcoG-2/CC2_data/filtered/SRR14295232_F_filt.fastq.gz

    ## Encountered 3082 unique sequences from 9584 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/EcoG-2/CC2_data/filtered/SRR14295233_F_filt.fastq.gz

    ## Encountered 2829 unique sequences from 7289 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/EcoG-2/CC2_data/filtered/SRR14295234_F_filt.fastq.gz

    ## Encountered 2262 unique sequences from 5733 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/EcoG-2/CC2_data/filtered/SRR14295235_F_filt.fastq.gz

    ## Encountered 3412 unique sequences from 10466 total sequences read.

``` r
derepRs <- derepFastq(filtRs, verbose=TRUE)
```

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/EcoG-2/CC2_data/filtered/SRR14295221_R_filt.fastq.gz

    ## Encountered 3594 unique sequences from 9498 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/EcoG-2/CC2_data/filtered/SRR14295222_R_filt.fastq.gz

    ## Encountered 1682 unique sequences from 5674 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/EcoG-2/CC2_data/filtered/SRR14295223_R_filt.fastq.gz

    ## Encountered 4019 unique sequences from 11310 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/EcoG-2/CC2_data/filtered/SRR14295224_R_filt.fastq.gz

    ## Encountered 2277 unique sequences from 7029 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/EcoG-2/CC2_data/filtered/SRR14295225_R_filt.fastq.gz

    ## Encountered 2391 unique sequences from 6435 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/EcoG-2/CC2_data/filtered/SRR14295226_R_filt.fastq.gz

    ## Encountered 2393 unique sequences from 7301 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/EcoG-2/CC2_data/filtered/SRR14295227_R_filt.fastq.gz

    ## Encountered 2333 unique sequences from 7010 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/EcoG-2/CC2_data/filtered/SRR14295228_R_filt.fastq.gz

    ## Encountered 2672 unique sequences from 9493 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/EcoG-2/CC2_data/filtered/SRR14295229_R_filt.fastq.gz

    ## Encountered 2459 unique sequences from 7002 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/EcoG-2/CC2_data/filtered/SRR14295230_R_filt.fastq.gz

    ## Encountered 2693 unique sequences from 8647 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/EcoG-2/CC2_data/filtered/SRR14295231_R_filt.fastq.gz

    ## Encountered 2290 unique sequences from 7932 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/EcoG-2/CC2_data/filtered/SRR14295232_R_filt.fastq.gz

    ## Encountered 2898 unique sequences from 9584 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/EcoG-2/CC2_data/filtered/SRR14295233_R_filt.fastq.gz

    ## Encountered 2590 unique sequences from 7289 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/EcoG-2/CC2_data/filtered/SRR14295234_R_filt.fastq.gz

    ## Encountered 2134 unique sequences from 5733 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/EcoG-2/CC2_data/filtered/SRR14295235_R_filt.fastq.gz

    ## Encountered 3151 unique sequences from 10466 total sequences read.

``` r
names(derepFs) <- sampleNames
names(derepRs) <- sampleNames
```

Estimation du taux d’erreur

``` r
errF <- learnErrors(filtFs, multithread=TRUE )
```

    ## 16856420 total bases in 120403 reads from 15 samples will be used for learning the error rates.

``` r
errR <- learnErrors(filtRs, multithread=TRUE)  
```

    ## 15652390 total bases in 120403 reads from 15 samples will be used for learning the error rates.

``` r
plotErrors(errF, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](CC2_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
plotErrors(errR, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](CC2_files/figure-gfm/unnamed-chunk-13-2.png)<!-- --> Les taux
d’erreur pour chaque transition possible (A→C, A→G, …) sont indiqués.
Les points sont les taux d’erreur observés pour chaque score de qualité
du consensus.  
\* La ligne noire montre les taux d’erreur estimés.  
\* La ligne rouge montre les taux d’erreur attendus selon la définition
nominale du Q-score.  
Ici, les taux d’erreur estimés (ligne noire) correspondent bien aux taux
observés (points), et les taux d’erreur diminuent avec l’augmentation de
la qualité, comme prévu. Tout semble raisonnable et nous poursuivons
avec confiance.

Inférence

``` r
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)  
```

    ## Sample 1 - 9498 reads in 3524 unique sequences.
    ## Sample 2 - 5674 reads in 1784 unique sequences.
    ## Sample 3 - 11310 reads in 3893 unique sequences.
    ## Sample 4 - 7029 reads in 2261 unique sequences.
    ## Sample 5 - 6435 reads in 2470 unique sequences.
    ## Sample 6 - 7301 reads in 2609 unique sequences.
    ## Sample 7 - 7010 reads in 2571 unique sequences.
    ## Sample 8 - 9493 reads in 2755 unique sequences.
    ## Sample 9 - 7002 reads in 2612 unique sequences.
    ## Sample 10 - 8647 reads in 2758 unique sequences.
    ## Sample 11 - 7932 reads in 2505 unique sequences.
    ## Sample 12 - 9584 reads in 3082 unique sequences.
    ## Sample 13 - 7289 reads in 2829 unique sequences.
    ## Sample 14 - 5733 reads in 2262 unique sequences.
    ## Sample 15 - 10466 reads in 3412 unique sequences.

``` r
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
```

    ## Sample 1 - 9498 reads in 3594 unique sequences.
    ## Sample 2 - 5674 reads in 1682 unique sequences.
    ## Sample 3 - 11310 reads in 4019 unique sequences.
    ## Sample 4 - 7029 reads in 2277 unique sequences.
    ## Sample 5 - 6435 reads in 2391 unique sequences.
    ## Sample 6 - 7301 reads in 2393 unique sequences.
    ## Sample 7 - 7010 reads in 2333 unique sequences.
    ## Sample 8 - 9493 reads in 2672 unique sequences.
    ## Sample 9 - 7002 reads in 2459 unique sequences.
    ## Sample 10 - 8647 reads in 2693 unique sequences.
    ## Sample 11 - 7932 reads in 2290 unique sequences.
    ## Sample 12 - 9584 reads in 2898 unique sequences.
    ## Sample 13 - 7289 reads in 2590 unique sequences.
    ## Sample 14 - 5733 reads in 2134 unique sequences.
    ## Sample 15 - 10466 reads in 3151 unique sequences.

Inspection du premier objet.

``` r
dadaFs[[1]] 
```

    ## dada-class: object describing DADA2 denoising results
    ## 355 sequence variants were inferred from 3524 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

# Construire le tableau des séquences et supprimer les chimères

Fusion des reads appariés.

``` r
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs)
```

Inspecter le data.frame de fusion du premier échantillon.

``` r
head(mergers[[1]])
```

    ##                                                                                                                                                                                                                                                        sequence
    ## 1 TACAGAGGTCCCAAGCGTTGTTCGGATTCACTGGGCGTAAAGGGTGCGTAGGCGGTCGGGTAAGTCTGACGTGAAATCTTCGAGCTCAACTCGGAAACTGCGTCGGATACTATTCGGCTCGAGGAATGGAGGGGAGACTGGAATACTTGGTGTAGCAGTGAAATGCGTAGATATCAAGTGGAACACCAGTGGCGAAGGCGAGTCTCTGGACATTTCCTGACGCTGAGGCACGAAAGCCAGGGGAGCAAACGGG
    ## 2 TACGGAGGGTGCAAGCGTTATCCGGATTCACTGGGTTTAAAGGGTGCGTAGGCGGGCAGGTAAGTCAGTGGTGAAATCTCCGGGCTTAACCCGGAAACTGCCGTTGATACTACTTGTCTTGAATATTGTGGAGGTAAGCGGAATATGTCATGTAGCGGTGAAATGCTTAGATATGACATAGAACACCAATTGCGAAGGCAGCTTACTACACAATGATTGACGCTGAGGCACGAAAGCGTGGGGATCAAACAGG
    ## 3 TACGGAGGGTGCAAGCGTTATCCGGATTCACTGGGTTTAAAGGGTGCGTAGGCGGGCAGGTAAGTCAGTGGTGAAATCTCCGAGCTTAACTCGGAAACTGCCGTTGATACTACTTGTCTTGAATATTGTGGAGGTAAGCGGAATATGTCATGTAGCGGTGAAATGCTTAGATATGACATAGAACACCAATTGCGAAGGCAGCTTACTACACAATGATTGACGCTGAGGCACGAAAGCGTGGGGATCAAACAGG
    ## 4 TACGGAGGGTGCAAGCGTTATCCGGATTCACTGGGTTTAAAGGGTGCGTAGGTGGGTATGTAAGTCAGTGGTGAAATCCCCGAGCTCAACTTGGGAACTGCCGTTGATACTATATATCTTGAATGTCGTAGAGGTAAGCGGAATATGTCATGTAGCGGTGAAATGCTTAGAGATGACATAGAACACCAATTGCGAAGGCAGCTTACTATGCGAATATTGACACTGAGGCACGAAAGCGTGGGTAGCAAACAGG
    ## 5 TACAGAGGTCCCAAGCGTTGTTCGGATTCACTGGGCGTAAAGGGTGCGTAGGCGGTCGGGTAAGTCTGACGTGAAATCTTCAAGCTCAACTTGGAAACTGCGTCGGATACTATTCGGCTCGAGGAATGGAGGGGAGACTGGAATACTTGGTGTAGCAGTGAAATGCGTAGATATCAAGTGGAACACCAGTGGCGAAGGCGAGTCTCTGGACATTTCCTGACGCTGAGGCACGAAAGCCAGGGGAGCAAACGGG
    ## 6 TACGTAGGGCGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGTGCGCAGGCGGATCTGTAAGTCAGGGGTGAAATCTCGCGGCTCAACCGCGAAACTGCCTTTGATACTGTGGATCTTGAGTTTGGGAGAGGTTGATGGAATTCCAGGTGTAGCGGTGAAATGCGTAGATATCTGGAAGAACACCAGTGGCGAAGGCGGTCAACTGGCCCATAACTGACGCTCATGCACGAAAGCGTGGGGATCAAACAGG
    ##   abundance forward reverse nmatch nmismatch nindel prefer accept
    ## 1       476       2       1     17         0      0      2   TRUE
    ## 2       386       3       2     17         0      0      2   TRUE
    ## 3       376       1       2     17         0      0      2   TRUE
    ## 4       341       5       3     17         0      0      1   TRUE
    ## 5       340       4       1     17         0      0      2   TRUE
    ## 6       246       6       4     17         0      0      1   TRUE

Construction d’une table d’ASV

``` r
seqtabAll <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])  
table(nchar(getSequences(seqtabAll)))  
```

    ## 
    ##  221  222  223  225  241  246  247  248  249  250  251  252  253  254  255  256 
    ##    2    2    2    1    2    2    3    2    2    4    4   68 1060   57    9    3 
    ##  257  258 
    ##    2    1

Ce tableau contient 1226 ASV, et les longueurs de nos séquences
fusionnées se situent toutes dans la plage attendue pour cet amplicon
V4, en majorité à 253 nucléotides.

Suppression des chimères.

``` r
seqtabNoC <- removeBimeraDenovo(seqtabAll, verbose=TRUE)
```

    ## Identified 102 bimeras out of 1226 input sequences.

``` r
1-ncol(seqtabNoC)/ncol(seqtabAll)
```

    ## [1] 0.08319739

102 chimères ont été supprimées, soit 8% de nos ASVs.

# Assigner la taxonomie

Acquisition du jeu d’entraînement Silva (*non-exécuté*)

``` bash
cd home/rstudio
wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz   
wget https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz
```

Assignement taxonomique

``` r
taxTab <- assignTaxonomy(seqtabNoC, "/home/rstudio/EcoG-2/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxTab <- addSpecies(taxTab, "/home/rstudio/EcoG-2/silva_species_assignment_v138.1.fa.gz")

unname(head(taxTab))
```

    ##      [,1]       [,2]             [,3]                  [,4]              
    ## [1,] "Bacteria" "Bacteroidota"   "Bacteroidia"         "Chitinophagales" 
    ## [2,] "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Xanthomonadales" 
    ## [3,] "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Enterobacterales"
    ## [4,] "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Xanthomonadales" 
    ## [5,] "Bacteria" "Bacteroidota"   "Bacteroidia"         "Chitinophagales" 
    ## [6,] "Bacteria" "Bacteroidota"   "Bacteroidia"         "Flavobacteriales"
    ##      [,5]                 [,6]               [,7]
    ## [1,] "Chitinophagaceae"   NA                 NA  
    ## [2,] "Rhodanobacteraceae" "Rhodanobacter"    NA  
    ## [3,] "Aeromonadaceae"     "Aeromonas"        NA  
    ## [4,] "Rhodanobacteraceae" "Dokdonella"       NA  
    ## [5,] "Chitinophagaceae"   NA                 NA  
    ## [6,] "Weeksellaceae"      "Chryseobacterium" NA

Arbre phylogénétique

``` r
seqs <- getSequences(seqtabNoC)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
        rearrangement = "stochastic", control = pml.control(trace = 0))
plot(fitGTR)
```

![](CC2_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

# Phyloseq

Classification des échantillons

``` r
#Après des heures de tentatives, je me suis résigné à assigner l'origine de chaque échantillon mannuellement...
sampleNumber <- as.integer(sapply(strsplit(sampleNames, "52"), `[`, 2)) 
sampleType <- data.frame(sampleNames)
sampleType$Types <- "Boues"
sampleType$Types[sampleNumber==28] <- "Déjections" 
sampleType$Types[sampleNumber==30] <- "Déjections" 
sampleType$Types[sampleNumber==31] <- "Déjections" 
sampleType$Types[sampleNumber==32] <- "Déjections" 
sampleType$Types[sampleNumber==35] <- "Déjections" 
sampleType$Types[sampleNumber==21] <- "Compost" 
sampleType$Types[sampleNumber==22] <- "Compost" 
sampleType$Types[sampleNumber==23] <- "Compost" 
sampleType$Types[sampleNumber==24] <- "Compost" 
sampleType$Types[sampleNumber==25] <- "Compost" 
rownames(sampleType) <- sampleNames
```

Combiner les données dans un objet phyloseq

``` r
ps <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=FALSE), 
               sample_data(sampleType), 
               tax_table(taxTab),phy_tree(fitGTR$tree))  
ps <- prune_samples(sample_names(ps) != "Mock", ps)  
ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 1124 taxa and 15 samples ]
    ## sample_data() Sample Data:       [ 15 samples by 2 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 1124 taxa by 7 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 1124 tips and 1122 internal nodes ]

On remplace les séquences par ASV.

``` r
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 1124 taxa and 15 samples ]
    ## sample_data() Sample Data:       [ 15 samples by 2 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 1124 taxa by 7 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 1124 tips and 1122 internal nodes ]
    ## refseq()      DNAStringSet:      [ 1124 reference sequences ]

# Courbes de raréfaction

``` r
rarecurve(seqtabNoC, step=10, ylab="ASVs", label=F)
```

![](CC2_files/figure-gfm/unnamed-chunk-26-1.png)<!-- --> Cela montre que
la profondeur de séquençage est optimale. Cependant, de nombreuses
tentatives n’ont pas suffi à colorer les courbes en fonction de
l’échantillon (boues, déjections, compost).

# Abondance différentielle

``` r
ps.rel = transform_sample_counts(ps, function(x) x/sum(x)*100)   # Agglomération des taxa
glom <- tax_glom(ps.rel, taxrank = 'Phylum', NArm = FALSE)
ps.melt <- psmelt(glom)
ps.melt$Phylum <- as.character(ps.melt$Phylum)

ps.melt <- ps.melt %>%
  group_by(Types, Phylum) %>%
  mutate(median=median(Abundance))  # Sélection des groupes avec median>1
keep <- unique(ps.melt$Phylum[ps.melt$median > 1])
ps.melt$Phylum[!(ps.melt$Phylum %in% keep)] <- "others"   #pour obtenir les mêmes lignes ensemble
ps.melt_sum <- ps.melt %>%
  group_by(Types,Phylum) %>%
  summarise(Abundance=sum(Abundance)/5)
```

    ## `summarise()` has grouped output by 'Types'. You can override using the `.groups` argument.

``` r
ps.melt_sum$Types <- factor(ps.melt_sum$Types, levels=c("Boues", "Déjections", "Compost")) #pour mettre les barres dans le bon ordre

ggplot(ps.melt_sum, aes(x = Types, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", aes(fill=Phylum)) + 
  labs(x="", y="%") +
  facet_wrap(~Types, scales= "free_x", nrow=1) +
  theme_classic() + 
  theme(strip.background = element_blank(), 
        axis.text.x.bottom = element_text(angle=45,vjust=1,hjust=1))
```

![](CC2_files/figure-gfm/unnamed-chunk-27-1.png)<!-- --> A la différence
de la publication d’origine, ici les Planctomycetota sont présents en
abondance supérieure à celle des Sumerlaeota.

# Alpha-diversité

``` r
#chao1(ps) 
plot_richness(ps, x="Types", measures="Chao1")
```

![](CC2_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

# Bêta-diversité

``` r
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
```

    ## Run 0 stress 0.00009407978 
    ## Run 1 stress 0.00009711485 
    ## ... Procrustes: rmse 0.1280181  max resid 0.1689617 
    ## Run 2 stress 0.00009123674 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1379452  max resid 0.1911845 
    ## Run 3 stress 0.0008693728 
    ## Run 4 stress 0.00008667921 
    ## ... New best solution
    ## ... Procrustes: rmse 0.02625075  max resid 0.03704369 
    ## Run 5 stress 0.3546489 
    ## Run 6 stress 0.00009322212 
    ## ... Procrustes: rmse 0.005084721  max resid 0.007433207 
    ## Run 7 stress 0.00009818141 
    ## ... Procrustes: rmse 0.06213774  max resid 0.08256831 
    ## Run 8 stress 0.330287 
    ## Run 9 stress 0.00009737697 
    ## ... Procrustes: rmse 0.01765278  max resid 0.02517663 
    ## Run 10 stress 0.00009430361 
    ## ... Procrustes: rmse 0.04194645  max resid 0.05750509 
    ## Run 11 stress 0.00008965387 
    ## ... Procrustes: rmse 0.09322377  max resid 0.1182402 
    ## Run 12 stress 0.001185915 
    ## Run 13 stress 0.00007975256 
    ## ... New best solution
    ## ... Procrustes: rmse 0.07561421  max resid 0.09922308 
    ## Run 14 stress 0.00009638623 
    ## ... Procrustes: rmse 0.02496443  max resid 0.03385009 
    ## Run 15 stress 0.00009381844 
    ## ... Procrustes: rmse 0.1126189  max resid 0.149425 
    ## Run 16 stress 0.00009396349 
    ## ... Procrustes: rmse 0.02274824  max resid 0.03080654 
    ## Run 17 stress 0.001322823 
    ## Run 18 stress 0.00009635608 
    ## ... Procrustes: rmse 0.09897282  max resid 0.1363372 
    ## Run 19 stress 0.00009425525 
    ## ... Procrustes: rmse 0.06251906  max resid 0.08660971 
    ## Run 20 stress 0.00009831433 
    ## ... Procrustes: rmse 0.01169614  max resid 0.01570907 
    ## *** No convergence -- monoMDS stopping criteria:
    ##      3: no. of iterations >= maxit
    ##     15: stress < smin
    ##      2: stress ratio > sratmax

    ## Warning in metaMDS(veganifyOTU(physeq), distance, ...): stress is (nearly) zero:
    ## you may have insufficient data

``` r
plot_ordination(ps.prop, ord.nmds.bray, color="Types", title="Bray-Curtis")
```

![](CC2_files/figure-gfm/unnamed-chunk-29-1.png)<!-- --> Les
échantillons sont bien regroupés par type de prélèvement et les trois
groupes sont très éloignés les uns des autres, comme on le retrouve sur
la figure de l’article. Il y a donc très peu d’ASVs en commun entre les
boues, les déjections, le compost. Après de nombreuses tentatives, je ne
suis pas parvenu à mettre la légende dans l’ordre (boues, déjections,
compost).

# Diagramme de Venn

``` r
venn(3, snames=c("Boues","Déjections","Compost"))
```

![](CC2_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

# Abondance des ASVs en commun
