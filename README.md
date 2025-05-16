Functional Alternative Splicing Events Score (FAScore)
================
chenli
2024/11/8

</br>

## What is FAScore?

FAScore is a machine learning framework for predicting functional
alternative splicing (AS) events within biological processes, such as
hematopoietic differentiation. It integrates dynamic expression patterns
of AS events and their host genes, along with sequence conservation and
structural details of the isoforms, to compute a functional score via
Random Forest (RF) modeling. Additionally, FAScore uses a Gaussian
Mixture Model (GMM) to classify AS events into functional categories.

The dynamic features used in FAScore include expression breadth (*R*),
expression specificity (*τ*), the expression relationship (Spearman’s
rank coefficient of correlation *ρ*) and *p* value (*P*<sub>*ρ*</sub>),
expression gradient (the linear fitting slope *β*) and *p* value
(*P*<sub>*β*</sub>). The dynamic score (DyScore) for genes or AS events
is defined as the mean of these scaled dynamic features, normalized to a
range from 0 to 1. Consequently, DyScore values range from -1 to 1,
where positive values indicate up-regulation, and negative values
indicate down-regulation. The structural features, sourced from the
APPRIS database, include cross-species sequence conservation score, the
functionally important residues score, domain integrity score,
trans-membrane helices score, signal peptide score, subcellular location
score, and the structural homologs and integrity score. The final
outputs of FAScore are the dynamic score (DyScore) and functional score
(FAScore) for all annotated AS events.

The package support the multiple species including human, mouse,
zebrafish, rat, pig, chimpanzee, chicken, cow, macaque, fruitfly, and
elegans. But if you only calculate the DyScore, you can use this package
in any species. The following figure shows the workflow of FAScore.

</br> 

<div align="center">
<img width = '900' height ='320' src ="https://github.com/LuChenLab/FAScore/blob/main/workflow.png"/>
</div>
<br>


## Contents

[Installation](#Installation)  
[Input](#Input)  
[Output](#Output)  
[Tutorial](#Tutorial)  
[SessionInfo](#SessionInfo)  
[Contact](#Contact)  
[Citation](#Citation)

## Installation

FAScore has been developed with `R 4.0.0` and the following packages are
needed to be installed in `R` (R scripts in FAScore automatically check
to see if an R library is already installed and then install those that
are needed. So no need for manual preinstallation!):

    randomForest (>= 4.0),
    SummarizedExperiment (>= 1.18),
    dplyr (>= 1.0),
    GenomicRanges (>= 1.0),
    GenomicFeatures (>= 1.0),
    ggplot2 (>= 3.0),
    ggrepel (>= 0.9),
    methods (>= 4.0),
    rtracklayer (>= 1.0),
    stringr (>= 1.0),
    tidyr (>= 1.0),
    IRanges (>= 2.0),
    stats (>= 4.0),
    reshape2 (>= 1.0),
    parallel (>= 4.0),
    S4Vectors (>= 0.26),
    cowplot (>= 1.0)
    mclust (>= 5.0)

To install FAScore, you have two options: either install directly from
GitHub or use the compressed source file:

``` r
# Install from GitHub if remotes package is not installed
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github("LuChenLab/FAScore/package")
```

Alternatively, you can install FAScore using the source file downloaded
from the repository :

``` bash
# Install FAScore from a downloaded source file
R CMD INSTALL FAScore_0.2.0.tar.gz
```

</br> </br>

## Input

#### PSI value of all AS events

A data frame with alternative splicing (AS) events or splice junctions
(SJs) as rows and sample IDs as columns. The percent spliced in (PSI)
values are calculated using tools like rMAST or ICAS, among others.

#### Gene expression

A data frame with gene IDs as rows and sample IDs as columns, containing
normalized gene expression counts. The columns should match those in the
PSI data frame.

#### Isoform expression

A data frame with transcript IDs as rows and sample IDs as columns,
containing normalized isoform expression counts (e.g., from RSEM). The
columns should be consistent with those in the PSI data frame.

#### Gene annotation file (GTF file)

The file path to a gene transfer format (GTF) file, which contains
genome annotation details for AS (SJs), isoforms, and gene loci. GTF
files specific to various species can be obtained from the
[**Ensembl**](https://asia.ensembl.org/info/data/ftp/index.html)
database.

#### Sample information

A data frame with sample IDs and group variables in the colData, where
rows correspond to columns in the PSI, gene expression, and isoform
expression data frames. The group variable must contain at least three
distinct classes.

#### Isoform structural information (optional)

A data frame with structural-related scores such as cross-species
sequence conservation score (Conservation), the functionally important
residues score (AA.residues), domain integrity score (Domain.integrity),
the structural homologs and integrity score (Protein3D), trans-membrane
helices score (Transmembrane), signal peptide score (SignalP), and
subcellular location score (TargetP). Files for species beyond
**human**(GRCh38), **rhesus macaque**(Mmul10), **rat**(Rnor6.0),
**mouse**(GRCm38), and **zebrafish**(GRCz11) can be downloaded from the
[**APPRIS**](https://appris.bioinfo.cnio.es/#/) database. </br> </br>

## Output

The output will be an S4 object containing the following slots:

-   `assays` : A SimpleList of matrix-like objects, including the PSI
    matrix of AS events derived from the input.
-   `colData` : An optional DataFrame providing information on the
    samples, as specified in the input.
-   `rowRanges` : A GRanges object describing the ranges of AS events
    and corresponding gene and isoform IDs. It includes columns for
    **ASID**, **TranscriptID**, **GeneID**, **GeneName** and
    **matchTransID**. The **matchTransID** column stores the
    best-matching transcript for each gene based on maximum expression
    values (from `ISOFORM` data), aiding in mapping isoform structural
    information.
-   `GENE` : A data frame of gene expression values provided in the
    input.
-   `ISOFORM` : A data frame of isoform expression values provided in
    the input.
-   `Correlation` : A list containing the Spearman’s rank correlation
    coefficient (*ρ*) and corresponding p-value (*P*<sub>*ρ*<sub>) for
    genes or AS events.
-   `Tau` : A list containing the specificity index (*τ*) for genes and
    AS events, which quantifies the specificity of expression profiles.
    Higher values indicate greater specificity.
-   `Range` : A list containing range value (*R*) for genes and AS
    events, representing the breadth of expression.
-   `Linear` : A list containing linear fit parameters (intercept, slope
    (*β*), r<sup>2</sup>), and the p-value (*P*<sub>*β*</sub>) for the
    linear fit of genes and AS events, representing the expression
    gradient.
-   `DyScore` : A list capturing dynamic scores for genes and AS events.
-   `Structure` : A data frame of isoform structural scores obtained
    from the APPRIS database.
-   `RFpredict` : The predicted scores and classification results from
    the pre-trained model, stored in the **FAScore** and **FAStype**
    columns.
-   `GMM` : The parameter result of gaussian mixture model (GMM).

</br> </br>

## Tutorial

### Usage

For detailed usage of parameters of functions, please `help()` or `?`.

### Prepare

Loading example data set.

``` r
library(FAScore)
library(dplyr)

data(ExDataSet)
```

``` r
head(ExDataSet$AS)
```

    ##                      LT-HSC_1   LT-HSC_2   ST-HSC_1   ST-HSC_2
    ## 1:4774517-4776409:-  1.485149   1.298701  7.1428571  5.9523810
    ## 1:4774517-4777524:- 98.514851  98.701299 92.8571429 94.0476190
    ## 1:4782734-4783950:- 98.007590  97.067449 97.0588235 99.3775934
    ## 1:4782734-4785572:-  1.992410   2.932551  2.9411765  0.6224066
    ## 1:4828650-4830267:+ 99.189463 100.000000 99.6891192 99.2125984
    ## 1:4828650-4832310:+  0.810537   0.000000  0.3108808  0.7874016
    ##                          MPP_1      MPP_2      CMP_1      CMP_2
    ## 1:4774517-4776409:-  3.1446541  5.1282051  3.7037037  5.8394161
    ## 1:4774517-4777524:- 96.8553459 94.8717949 96.2962963 94.1605839
    ## 1:4782734-4783950:- 98.2180294 97.3895582 98.0341880 96.0591133
    ## 1:4782734-4785572:-  1.7819706  2.6104418  1.9658120  3.9408867
    ## 1:4828650-4830267:+ 99.7979798 99.4094488 99.4541485 99.0825688
    ## 1:4828650-4832310:+  0.2020202  0.5905512  0.5458515  0.9174312
    ##                          MEP_1     MEP_2       MK_1       MK_2
    ## 1:4774517-4776409:-  3.3519553  4.901961  2.0408163 10.2564103
    ## 1:4774517-4777524:- 96.6480447 95.098039 97.9591837 89.7435897
    ## 1:4782734-4783950:- 98.4264786 97.285714 98.1351981 98.7500000
    ## 1:4782734-4785572:-  1.5735214  2.714286  1.8648019  1.2500000
    ## 1:4828650-4830267:+ 99.1423671 98.679245 99.3036212 99.7732426
    ## 1:4828650-4832310:+  0.8576329  1.320755  0.6963788  0.2267574
    ##                         EryA_1     EryA_2     EryB_1     EryB_2
    ## 1:4774517-4776409:-  3.4090909  1.6393443  13.333333  16.666667
    ## 1:4774517-4777524:- 96.5909091 98.3606557  86.666667  83.333333
    ## 1:4782734-4783950:- 97.8494624 99.6721311  98.648649  98.461538
    ## 1:4782734-4785572:-  2.1505376  0.3278689   1.351351   1.538462
    ## 1:4828650-4830267:+ 99.0291262 98.9510490 100.000000 100.000000
    ## 1:4828650-4832310:+  0.9708738  1.0489510   0.000000   0.000000

``` r
head(ExDataSet$GENE)
```

    ##                      LT-HSC_1 LT-HSC_2 ST-HSC_1 ST-HSC_2 MPP_1 MPP_2
    ## ENSMUSG00000033845      23.05    13.37    24.63    16.83 18.55 17.75
    ## ENSMUSG00000025903      25.01    16.63    27.51    20.29 23.53 20.05
    ## ENSMUSG00000104217       0.18     0.04     0.29     0.04  0.43  0.15
    ## ENSMUSG00000025903.1    25.01    16.63    27.51    20.29 23.53 20.05
    ## ENSMUSG00000104217.1     0.18     0.04     0.29     0.04  0.43  0.15
    ## ENSMUSG00000033813      40.78    30.87    43.12    33.82 34.04 32.50
    ##                      CMP_1 CMP_2 MEP_1 MEP_2  MK_1  MK_2 EryA_1
    ## ENSMUSG00000033845   17.00 14.88 31.98 24.67  7.66 12.33   9.95
    ## ENSMUSG00000025903   21.56 14.42 13.68  8.56 14.67 13.55   7.89
    ## ENSMUSG00000104217    0.15  0.04  0.07  0.14  0.05  0.15   0.06
    ## ENSMUSG00000025903.1 21.56 14.42 13.68  8.56 14.67 13.55   7.89
    ## ENSMUSG00000104217.1  0.15  0.04  0.07  0.14  0.05  0.15   0.06
    ## ENSMUSG00000033813   37.95 20.42 45.18 27.30 12.08 17.55  32.63
    ##                      EryA_2 EryB_1 EryB_2
    ## ENSMUSG00000033845    10.21   2.18   1.65
    ## ENSMUSG00000025903     9.41   1.58   1.52
    ## ENSMUSG00000104217     0.22   0.03   0.00
    ## ENSMUSG00000025903.1   9.41   1.58   1.52
    ## ENSMUSG00000104217.1   0.22   0.03   0.00
    ## ENSMUSG00000033813    45.86  24.17  21.21

``` r
head(ExDataSet$ISOFORM)
```

    ##                               gene_id LT-HSC_1 LT-HSC_2 ST-HSC_1
    ## ENSMUST00000001166 ENSMUSG00000001138     5.13     4.83     5.08
    ## ENSMUST00000001171 ENSMUSG00000001143     0.68     2.06     1.78
    ## ENSMUST00000001172 ENSMUSG00000079610     1.07     0.43     1.00
    ## ENSMUST00000003219 ENSMUSG00000003135     8.00     6.05     8.05
    ## ENSMUST00000006037 ENSMUSG00000005886     6.71     5.55     6.82
    ## ENSMUST00000006462 ENSMUSG00000006299    32.24    22.48    36.18
    ##                    ST-HSC_2 MPP_1 MPP_2 CMP_1 CMP_2 MEP_1 MEP_2
    ## ENSMUST00000001166     5.49  4.10  4.73  3.71  2.80  2.63  2.19
    ## ENSMUST00000001171     1.52  0.97  1.46  0.56  0.78  0.21  0.26
    ## ENSMUST00000001172     0.69  1.27  0.62  0.76  0.73  0.69  0.44
    ## ENSMUST00000003219     6.85  6.38  6.81  5.82  3.67  5.99  3.72
    ## ENSMUST00000006037     5.47  5.55  4.64  5.13  3.25  3.59  2.65
    ## ENSMUST00000006462    28.38 29.68 26.82 24.80 18.09 25.18 19.33
    ##                     MK_1  MK_2 EryA_1 EryA_2 EryB_1 EryB_2
    ## ENSMUST00000001166  2.47  4.38   1.56   2.27   1.03   0.76
    ## ENSMUST00000001171  1.45  1.47   0.68   1.07   0.23   0.11
    ## ENSMUST00000001172  0.37  0.77   0.34   0.27   0.04   0.05
    ## ENSMUST00000003219  2.56  4.52   2.58   3.71   0.89   0.79
    ## ENSMUST00000006037  3.17  2.23   1.48   2.23   0.90   0.98
    ## ENSMUST00000006462 13.12 22.40   8.30  10.13   2.67   2.33

``` r
head(ExDataSet$meta)
```

    ##               Run CellType    Lineage
    ## LT-HSC_1 LT-HSC_1   LT-HSC       Stem
    ## LT-HSC_2 LT-HSC_2   LT-HSC       Stem
    ## ST-HSC_1 ST-HSC_1   ST-HSC       Stem
    ## ST-HSC_2 ST-HSC_2   ST-HSC       Stem
    ## MPP_1       MPP_1      MPP Progenitor
    ## MPP_2       MPP_2      MPP Progenitor

</br>

### Creating a FAScore object

The alternative splicing PSI matrix, gene expression matrix, transcripts
expression matrix, and sample information table need to be prepared to
create a FAScore object.

``` r
MyObj <- FAScoreDataSet(colData = ExDataSet$meta, AS = ExDataSet$AS, 
                        GENE = ExDataSet$GENE, ISOFORM = ExDataSet$ISOFORM)
```

``` r
MyObj
```

    ## class: FAScore 
    ## dim: 1000 16 
    ## metadata(1): version
    ## assays(1): AS
    ## rownames(1000): 1:4774517-4776409:- 1:4774517-4777524:- ...
    ##   1:86583097-86587015:+ 1:86583097-86589163:+
    ## rowData names(0):
    ## colnames(16): LT-HSC_1 LT-HSC_2 ... EryB_1 EryB_2
    ## colData names(3): Run CellType Lineage

</br>

### Calculate the dynamic scores of AS events and genes

The dynamic score of genes or AS events is defined as the mean of the
scaled values (ranging from 0 to 1) of dynamic features calculated by
the `CalcuFeature` function. These features include R<sub>scaled</sub>,
*τ*, *ρ*, *P*<sub>*ρ*</sub>, *β* and *P*<sub>*β*</sub>.
R<sub>scaled</sub> represents the range values of TPM and PSI across all
stages or cell types per gene or AS event, scaled to a range between 0
and 1 by defining the maximum range (max range). The specificity index
*τ* quantifies the specificity of an expression profile for a particular
stage or cell type, ranging from 0 (housekeeping) to 1 (stage or cell
type-specific). *ρ* is the Spearman’s rank correlation coefficient, and
*P*<sub>*ρ*</sub> is the corresponding *P*-value (two. side t-test). If
the *P*<sub>*ρ*</sub> is significant, it is set to 1, otherwise to 0.
*β* is the slope of the linear fit based on TPM and PSI across all
stages or cell types per gene or AS event, scaled between -1 and 1.
Values greater than or equal to 1 are set to 1, and values less than or
equal to -1 are set to -1. *P*<sub>*β*</sub> is the p-value from the
F-statistic test during the linear fitting; if the p-value is
significant, it is set to 1, otherwise to 0. The absolute DyScore ranges
between 0 and 1, with the sign indicating the direction of the slope
(*β*). The dynamic score is calculated using the `CalcuDyScore`
function.

</br>

Preprocessing the related dynamic features

``` r
MyObj <- CalcuFeature(MyObj, group.by = "CellType", cores = 10)
```

``` r
head(MyObj@Correlation$Gene)  # head(MyObj@Correlation$AS)
```

    ##                      Spearman.cor   Spearman.p
    ## ENSMUSG00000033845     -0.6449817 0.0069804230
    ## ENSMUSG00000025903     -0.8698377 0.0000119491
    ## ENSMUSG00000104217     -0.3452381 0.1903171729
    ## ENSMUSG00000025903.1   -0.8698377 0.0000119491
    ## ENSMUSG00000104217.1   -0.3452381 0.1903171729
    ## ENSMUSG00000033813     -0.3372840 0.2014142966

``` r
head(MyObj@Tau$Gene)  # head(MyObj@Tau$AS)
```

    ##   ENSMUSG00000033845   ENSMUSG00000025903   ENSMUSG00000104217 
    ##            0.5207666            0.4260012            0.6403941 
    ## ENSMUSG00000025903.1 ENSMUSG00000104217.1   ENSMUSG00000033813 
    ##            0.4260012            0.6403941            0.2337695

``` r
head(MyObj@Range$Gene)  # head(MyObj@Range$AS)
```

    ##   ENSMUSG00000033845   ENSMUSG00000025903   ENSMUSG00000104217 
    ##               26.410               22.350                0.275 
    ## ENSMUSG00000025903.1 ENSMUSG00000104217.1   ENSMUSG00000033813 
    ##               22.350                0.275               24.430

``` r
head(MyObj@Linear$Gene)  # head(MyObj@Linear$AS)
```

    ##                      intercept       slope        r2       pvalue
    ## ENSMUSG00000033845    42.13620 -0.77794989 0.3170984 0.0231327141
    ## ENSMUSG00000025903    44.48707 -1.09057858 0.5728384 0.0006882753
    ## ENSMUSG00000104217     2.13947 -0.06570322 0.1216056 0.1855878360
    ## ENSMUSG00000025903.1  44.48707 -1.09057858 0.5728384 0.0006882753
    ## ENSMUSG00000104217.1   2.13947 -0.06570322 0.1216056 0.1855878360
    ## ENSMUSG00000033813    48.20751 -0.26908852 0.1446797 0.1461264931

</br>

Calculating the dynamic scores

``` r
MyObj <- CalcuDyScore(MyObj, maxRange = 10, maxSlope = 1, type = "Gene")
```

``` r
MyObj <- CalcuDyScore(MyObj, maxRange = 30, maxSlope = 1, type = "AS")
```

``` r
head(MyObj@DyScore$Gene)
```

    ##   ENSMUSG00000033845   ENSMUSG00000025903   ENSMUSG00000104217 
    ##           -0.8239497           -0.8826398           -0.1752226 
    ## ENSMUSG00000025903.1 ENSMUSG00000104217.1   ENSMUSG00000033813 
    ##           -0.8826398           -0.1752226           -0.3066903

``` r
head(MyObj@DyScore$AS)
```

    ## 1:4774517-4776409:- 1:4774517-4777524:- 1:4782734-4783950:- 
    ##          0.29129371         -0.17977795          0.08980269 
    ## 1:4782734-4785572:- 1:4828650-4830267:+ 1:4828650-4832310:+ 
    ##         -0.15165456         -0.01005425          0.09792982

</br>

### Calculate the functional scores of AS events

We use dynamic, conserved and structural features to train a predictive
model using the `randomForest` function (pre-trained model). The
functional scores and classes can be calculated using the `CalcuFAScore`
function, which applies either the pre-trained model or a custom model
of your choice.

</br>

Adding gene and transcript information for AS events. Note that
processing the `HostGene` may take a considerable amount of time if the
GTF file is large.

``` r
my_gtf <- system.file("extdata", "GRCm38.93_sub.gtf.gz", package = "FAScore", mustWork = TRUE)

MyObj <- ASmapIso(MyObj, gtf = my_gtf, AStype = "exonic", cores = 10)
MyObj <- ChooseIso(MyObj)
```

</br>

Adding the structural scores of transcripts

``` r
MyObj <- matchAppris(MyObj, gtf = my_gtf, species = "MusMus")

head(MyObj@Structure)
```

</br>

or other species and annotation versions downloaded from
[**APPRIS**](https://appris.bioinfo.cnio.es/#/) database

``` r
MyObj <- matchAppris(MyObj, gtf = my_gtf, appris = OtherFile)

head(MyObj@Structure)
```

</br>

Calculating the functional scores and classes

``` r
MyObj <- CalcuFAScore(MyObj)
```

</br>

Sorted by FAScore

``` r
MyObj@RFpredict <- MyObj@RFpredict[order(MyObj@RFpredict$FAScore,decreasing = T),]

MyObj@RFpredict[,c(1,3,4,26,27)] %>% head
```

</br>

### Visualization

The `plotDyScore` function can rank the dynamic scores of genes and AS
events. You can label different classes of AS events, as defined by you,
along with their corresponding genes in the figure.

``` r
ASID = list(Func = (MyObj@RFpredict$AS %>% head), nonFunc  = (MyObj@RFpredict$AS %>% tail))

label = paste0(c(head(MyObj@RFpredict$GeneName),tail(MyObj@RFpredict$GeneName)), "|", stringr::str_split_fixed(unlist(ASID), "[|]", 2)[,2] )
```

``` r
plotDyScore(MyObj, ASID = ASID, label = label, color = c( "#FF4500", "#2E9FDF"), size = 18)
```

</br>

The `plotFAScore` function can be used to rank the functional scores of
AS events, with the option to label specific events of interest.

``` r
plotFAScore(MyObj, ASID = ASID, label = label, color = c( "#FF4500", "#2E9FDF"), size = 18) 
```

</br>  
</br>

## SessionInfo

``` r
sessionInfo()
```

    ## R version 4.1.2 (2021-11-01)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.6 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/atlas/libblas.so.3.10.3
    ## LAPACK: /usr/lib/x86_64-linux-gnu/atlas/liblapack.so.3.10.3
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods  
    ## [7] base     
    ## 
    ## other attached packages:
    ## [1] dplyr_1.1.4   FAScore_0.2.0 rmarkdown_2.9
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] readxl_1.3.1                shadowtext_0.0.8           
    ##   [3] backports_1.2.1             Hmisc_4.5-0                
    ##   [5] fastmatch_1.1-3             BiocFileCache_1.99.3       
    ##   [7] plyr_1.8.6                  igraph_1.2.6               
    ##   [9] lazyeval_0.2.2              splines_4.1.2              
    ##  [11] BiocParallel_1.25.5         GenomeInfoDb_1.27.11       
    ##  [13] ggplot2_3.5.1               digest_0.6.27              
    ##  [15] yulab.utils_0.1.1           htmltools_0.5.1.1          
    ##  [17] GOSemSim_2.17.1             viridis_0.6.1              
    ##  [19] GO.db_3.12.1                fansi_0.5.0                
    ##  [21] magrittr_2.0.1              checkmate_2.0.0            
    ##  [23] memoise_2.0.0               cluster_2.1.2              
    ##  [25] openxlsx_4.2.4              Biostrings_2.59.2          
    ##  [27] graphlayouts_0.7.1          matrixStats_0.62.0         
    ##  [29] prettyunits_1.1.1           enrichplot_1.11.2          
    ##  [31] jpeg_0.1-9                  colorspace_2.0-2           
    ##  [33] rappdirs_0.3.3              blob_1.2.2                 
    ##  [35] ggrepel_0.9.4               haven_2.4.1                
    ##  [37] xfun_0.24                   RCurl_1.98-1.3             
    ##  [39] crayon_1.4.1                jsonlite_1.7.2             
    ##  [41] scatterpie_0.1.6            survival_3.2-13            
    ##  [43] ape_5.5                     glue_1.6.2                 
    ##  [45] polyclip_1.10-0             gtable_0.3.0               
    ##  [47] zlibbioc_1.37.0             XVector_0.31.1             
    ##  [49] DelayedArray_0.17.10        car_3.0-11                 
    ##  [51] BiocGenerics_0.37.1         abind_1.4-5                
    ##  [53] scales_1.3.0                DOSE_3.17.0                
    ##  [55] DBI_1.1.1                   rstatix_0.7.0              
    ##  [57] Rcpp_1.0.12                 progress_1.2.2             
    ##  [59] viridisLite_0.4.0           htmlTable_2.2.1            
    ##  [61] tidytree_0.4.6              mclust_5.4.9               
    ##  [63] foreign_0.8-81              bit_4.0.4                  
    ##  [65] Formula_1.2-4               stats4_4.1.2               
    ##  [67] htmlwidgets_1.5.3           httr_1.4.2                 
    ##  [69] fgsea_1.17.0                RColorBrewer_1.1-2         
    ##  [71] ellipsis_0.3.2              XML_3.99-0.6               
    ##  [73] pkgconfig_2.0.3             farver_2.1.0               
    ##  [75] dbplyr_2.1.1                nnet_7.3-16                
    ##  [77] utf8_1.2.2                  tidyselect_1.2.0           
    ##  [79] rlang_1.1.2                 reshape2_1.4.4             
    ##  [81] AnnotationDbi_1.53.1        munsell_0.5.0              
    ##  [83] cellranger_1.1.0            tools_4.1.2                
    ##  [85] cachem_1.0.5                downloader_0.4             
    ##  [87] cli_3.6.1                   generics_0.1.0             
    ##  [89] RSQLite_2.2.7               broom_0.7.9                
    ##  [91] evaluate_0.14               stringr_1.5.1              
    ##  [93] fastmap_1.1.0               yaml_2.2.1                 
    ##  [95] ggtree_3.10.0               knitr_1.33                 
    ##  [97] bit64_4.0.5                 fs_1.5.0                   
    ##  [99] tidygraph_1.2.0             zip_2.2.0                  
    ## [101] randomForest_4.6-14         purrr_1.0.2                
    ## [103] KEGGREST_1.31.1             ggraph_2.0.5               
    ## [105] nlme_3.1-152                aplot_0.0.6                
    ## [107] DO.db_2.9                   biomaRt_2.47.7             
    ## [109] compiler_4.1.2              rstudioapi_0.13            
    ## [111] filelock_1.0.2              curl_4.3.2                 
    ## [113] png_0.1-7                   ggsignif_0.6.4             
    ## [115] treeio_1.15.6               tibble_3.2.1               
    ## [117] tweenr_1.0.2                stringi_1.7.3              
    ## [119] GenomicFeatures_1.43.8      forcats_0.5.1              
    ## [121] lattice_0.20-45             Matrix_1.6-0               
    ## [123] vctrs_0.6.4                 pillar_1.9.0               
    ## [125] lifecycle_1.0.3             BiocManager_1.30.22        
    ## [127] bitops_1.0-7                data.table_1.14.0          
    ## [129] cowplot_1.1.1               GenomicRanges_1.43.4       
    ## [131] rtracklayer_1.51.5          patchwork_1.1.1            
    ## [133] qvalue_2.23.0               BiocIO_1.1.2               
    ## [135] R6_2.5.0                    latticeExtra_0.6-29        
    ## [137] gridExtra_2.3               rio_0.5.27                 
    ## [139] IRanges_2.25.7              assertthat_0.2.1           
    ## [141] MASS_7.3-54                 SummarizedExperiment_1.21.3
    ## [143] rjson_0.2.20                GenomicAlignments_1.27.2   
    ## [145] Rsamtools_2.7.2             GenomeInfoDbData_1.2.4     
    ## [147] S4Vectors_0.29.15           parallel_4.1.2             
    ## [149] hms_1.1.0                   clusterProfiler_3.99.1     
    ## [151] grid_4.1.2                  rpart_4.1-15               
    ## [153] ggfun_0.1.3                 tidyr_1.3.0                
    ## [155] rvcheck_0.1.8               MatrixGenerics_1.3.1       
    ## [157] carData_3.0-4               ggpubr_0.4.0               
    ## [159] ggforce_0.3.3               Biobase_2.51.0             
    ## [161] base64enc_0.1-3             restfulr_0.0.13

</br>  
</br> </br>

## Contact

Please contact Lu Chen (<luchen@scu.edu.cn>) or Li Chen
(<chenli5679@163.com>).

</br>

## Citation

If you use FAScore in your publication, please cite FAScore by
