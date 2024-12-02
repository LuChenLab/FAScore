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

</br> </br>

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

    # Install from GitHub if remotes package is not installed
    if (!requireNamespace("remotes", quietly = TRUE))
        install.packages("remotes")

    remotes::install_github("LuChenLab/DEMINERS/DecodeR")

Alternatively, you can install FAScore using the source file downloaded
from the repository :

    # Install DecodeR from a downloaded source file
    R CMD INSTALL FAScore_0.2.0.tar.gz

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

    library(FAScore)
    library(dplyr)

    data(ExDataSet)

    head(ExDataSet$AS)

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

    head(ExDataSet$GENE)

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

    head(ExDataSet$ISOFORM)

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

    head(ExDataSet$meta)

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

    MyObj <- FAScoreDataSet(colData = ExDataSet$meta, AS = ExDataSet$AS, 
                            GENE = ExDataSet$GENE, ISOFORM = ExDataSet$ISOFORM)

    MyObj

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

    MyObj <- CalcuFeature(MyObj, group.by = "CellType", cores = 10)

    head(MyObj@Correlation$Gene)  # head(MyObj@Correlation$AS)

    ##                      Spearman.cor   Spearman.p
    ## ENSMUSG00000033845     -0.6449817 0.0069804230
    ## ENSMUSG00000025903     -0.8698377 0.0000119491
    ## ENSMUSG00000104217     -0.3452381 0.1903171729
    ## ENSMUSG00000025903.1   -0.8698377 0.0000119491
    ## ENSMUSG00000104217.1   -0.3452381 0.1903171729
    ## ENSMUSG00000033813     -0.3372840 0.2014142966

    head(MyObj@Tau$Gene)  # head(MyObj@Tau$AS)

    ##   ENSMUSG00000033845   ENSMUSG00000025903   ENSMUSG00000104217 
    ##            0.5207666            0.4260012            0.6403941 
    ## ENSMUSG00000025903.1 ENSMUSG00000104217.1   ENSMUSG00000033813 
    ##            0.4260012            0.6403941            0.2337695

    head(MyObj@Range$Gene)  # head(MyObj@Range$AS)

    ##   ENSMUSG00000033845   ENSMUSG00000025903   ENSMUSG00000104217 
    ##               26.410               22.350                0.275 
    ## ENSMUSG00000025903.1 ENSMUSG00000104217.1   ENSMUSG00000033813 
    ##               22.350                0.275               24.430

    head(MyObj@Linear$Gene)  # head(MyObj@Linear$AS)

    ##                      intercept       slope        r2       pvalue
    ## ENSMUSG00000033845    42.13620 -0.77794989 0.3170984 0.0231327141
    ## ENSMUSG00000025903    44.48707 -1.09057858 0.5728384 0.0006882753
    ## ENSMUSG00000104217     2.13947 -0.06570322 0.1216056 0.1855878360
    ## ENSMUSG00000025903.1  44.48707 -1.09057858 0.5728384 0.0006882753
    ## ENSMUSG00000104217.1   2.13947 -0.06570322 0.1216056 0.1855878360
    ## ENSMUSG00000033813    48.20751 -0.26908852 0.1446797 0.1461264931

</br>

Calculating the dynamic scores

    MyObj <- CalcuDyScore(MyObj, maxRange = 10, maxSlope = 1, type = "Gene")

    ## [1] 1
    ## [1] 2
    ## [1] 3
    ## [1] 4
    ## [1] 5
    ## [1] 6
    ## [1] 7
    ## [1] 8
    ## [1] 9
    ## [1] 10
    ## [1] 11
    ## [1] 12
    ## [1] 13
    ## [1] 14
    ## [1] 15
    ## [1] 16
    ## [1] 17
    ## [1] 18
    ## [1] 19
    ## [1] 20
    ## [1] 21
    ## [1] 22
    ## [1] 23
    ## [1] 24
    ## [1] 25
    ## [1] 26
    ## [1] 27
    ## [1] 28
    ## [1] 29
    ## [1] 30
    ## [1] 31
    ## [1] 32
    ## [1] 33
    ## [1] 34
    ## [1] 35
    ## [1] 36
    ## [1] 37
    ## [1] 38
    ## [1] 39
    ## [1] 40
    ## [1] 41
    ## [1] 42
    ## [1] 43
    ## [1] 44
    ## [1] 45
    ## [1] 46
    ## [1] 47
    ## [1] 48
    ## [1] 49
    ## [1] 50
    ## [1] 51
    ## [1] 52
    ## [1] 53
    ## [1] 54
    ## [1] 55
    ## [1] 56
    ## [1] 57
    ## [1] 58
    ## [1] 59
    ## [1] 60
    ## [1] 61
    ## [1] 62
    ## [1] 63
    ## [1] 64
    ## [1] 65
    ## [1] 66
    ## [1] 67
    ## [1] 68
    ## [1] 69
    ## [1] 70
    ## [1] 71
    ## [1] 72
    ## [1] 73
    ## [1] 74
    ## [1] 75
    ## [1] 76
    ## [1] 77
    ## [1] 78
    ## [1] 79
    ## [1] 80
    ## [1] 81
    ## [1] 82
    ## [1] 83
    ## [1] 84
    ## [1] 85
    ## [1] 86
    ## [1] 87
    ## [1] 88
    ## [1] 89
    ## [1] 90
    ## [1] 91
    ## [1] 92
    ## [1] 93
    ## [1] 94
    ## [1] 95
    ## [1] 96
    ## [1] 97
    ## [1] 98
    ## [1] 99
    ## [1] 100
    ## [1] 101
    ## [1] 102
    ## [1] 103
    ## [1] 104
    ## [1] 105
    ## [1] 106
    ## [1] 107
    ## [1] 108
    ## [1] 109
    ## [1] 110
    ## [1] 111
    ## [1] 112
    ## [1] 113
    ## [1] 114
    ## [1] 115
    ## [1] 116
    ## [1] 117
    ## [1] 118
    ## [1] 119
    ## [1] 120
    ## [1] 121
    ## [1] 122
    ## [1] 123
    ## [1] 124
    ## [1] 125
    ## [1] 126
    ## [1] 127
    ## [1] 128
    ## [1] 129
    ## [1] 130
    ## [1] 131
    ## [1] 132
    ## [1] 133
    ## [1] 134
    ## [1] 135
    ## [1] 136
    ## [1] 137
    ## [1] 138
    ## [1] 139
    ## [1] 140
    ## [1] 141
    ## [1] 142
    ## [1] 143
    ## [1] 144
    ## [1] 145
    ## [1] 146
    ## [1] 147
    ## [1] 148
    ## [1] 149
    ## [1] 150
    ## [1] 151
    ## [1] 152
    ## [1] 153
    ## [1] 154
    ## [1] 155
    ## [1] 156
    ## [1] 157
    ## [1] 158
    ## [1] 159
    ## [1] 160
    ## [1] 161
    ## [1] 162
    ## [1] 163
    ## [1] 164
    ## [1] 165
    ## [1] 166
    ## [1] 167
    ## [1] 168
    ## [1] 169
    ## [1] 170
    ## [1] 171
    ## [1] 172
    ## [1] 173
    ## [1] 174
    ## [1] 175
    ## [1] 176
    ## [1] 177
    ## [1] 178
    ## [1] 179
    ## [1] 180
    ## [1] 181
    ## [1] 182
    ## [1] 183
    ## [1] 184
    ## [1] 185
    ## [1] 186
    ## [1] 187
    ## [1] 188
    ## [1] 189
    ## [1] 190
    ## [1] 191
    ## [1] 192
    ## [1] 193
    ## [1] 194
    ## [1] 195
    ## [1] 196
    ## [1] 197
    ## [1] 198
    ## [1] 199

    MyObj <- CalcuDyScore(MyObj, maxRange = 30, maxSlope = 1, type = "AS")

    ## [1] 1
    ## [1] 2
    ## [1] 3
    ## [1] 4
    ## [1] 5
    ## [1] 6
    ## [1] 7
    ## [1] 8
    ## [1] 9
    ## [1] 10
    ## [1] 11
    ## [1] 12
    ## [1] 13
    ## [1] 14
    ## [1] 15
    ## [1] 16
    ## [1] 17
    ## [1] 18
    ## [1] 19
    ## [1] 20
    ## [1] 21
    ## [1] 22
    ## [1] 23
    ## [1] 24
    ## [1] 25
    ## [1] 26
    ## [1] 27
    ## [1] 28
    ## [1] 29
    ## [1] 30
    ## [1] 31
    ## [1] 32
    ## [1] 33
    ## [1] 34
    ## [1] 35
    ## [1] 36
    ## [1] 37
    ## [1] 38
    ## [1] 39
    ## [1] 40
    ## [1] 41
    ## [1] 42
    ## [1] 43
    ## [1] 44
    ## [1] 45
    ## [1] 46
    ## [1] 47
    ## [1] 48
    ## [1] 49
    ## [1] 50
    ## [1] 51
    ## [1] 52
    ## [1] 53
    ## [1] 54
    ## [1] 55
    ## [1] 56
    ## [1] 57
    ## [1] 58
    ## [1] 59
    ## [1] 60
    ## [1] 61
    ## [1] 62
    ## [1] 63
    ## [1] 64
    ## [1] 65
    ## [1] 66
    ## [1] 67
    ## [1] 68
    ## [1] 69
    ## [1] 70
    ## [1] 71
    ## [1] 72
    ## [1] 73
    ## [1] 74
    ## [1] 75
    ## [1] 76
    ## [1] 77
    ## [1] 78
    ## [1] 79
    ## [1] 80
    ## [1] 81
    ## [1] 82
    ## [1] 83
    ## [1] 84
    ## [1] 85
    ## [1] 86
    ## [1] 87
    ## [1] 88
    ## [1] 89
    ## [1] 90
    ## [1] 91
    ## [1] 92
    ## [1] 93
    ## [1] 94
    ## [1] 95
    ## [1] 96
    ## [1] 97
    ## [1] 98
    ## [1] 99
    ## [1] 100
    ## [1] 101
    ## [1] 102
    ## [1] 103
    ## [1] 104
    ## [1] 105
    ## [1] 106
    ## [1] 107
    ## [1] 108
    ## [1] 109
    ## [1] 110
    ## [1] 111
    ## [1] 112
    ## [1] 113
    ## [1] 114
    ## [1] 115
    ## [1] 116
    ## [1] 117
    ## [1] 118
    ## [1] 119
    ## [1] 120
    ## [1] 121
    ## [1] 122
    ## [1] 123
    ## [1] 124
    ## [1] 125
    ## [1] 126
    ## [1] 127
    ## [1] 128
    ## [1] 129
    ## [1] 130
    ## [1] 131
    ## [1] 132
    ## [1] 133
    ## [1] 134
    ## [1] 135
    ## [1] 136
    ## [1] 137
    ## [1] 138
    ## [1] 139
    ## [1] 140
    ## [1] 141
    ## [1] 142
    ## [1] 143
    ## [1] 144
    ## [1] 145
    ## [1] 146
    ## [1] 147
    ## [1] 148
    ## [1] 149
    ## [1] 150
    ## [1] 151
    ## [1] 152
    ## [1] 153
    ## [1] 154
    ## [1] 155
    ## [1] 156
    ## [1] 157
    ## [1] 158
    ## [1] 159
    ## [1] 160
    ## [1] 161
    ## [1] 162
    ## [1] 163
    ## [1] 164
    ## [1] 165
    ## [1] 166
    ## [1] 167
    ## [1] 168
    ## [1] 169
    ## [1] 170
    ## [1] 171
    ## [1] 172
    ## [1] 173
    ## [1] 174
    ## [1] 175
    ## [1] 176
    ## [1] 177
    ## [1] 178
    ## [1] 179
    ## [1] 180
    ## [1] 181
    ## [1] 182
    ## [1] 183
    ## [1] 184
    ## [1] 185
    ## [1] 186
    ## [1] 187
    ## [1] 188
    ## [1] 189
    ## [1] 190
    ## [1] 191
    ## [1] 192
    ## [1] 193
    ## [1] 194
    ## [1] 195
    ## [1] 196
    ## [1] 197
    ## [1] 198
    ## [1] 199
    ## [1] 200
    ## [1] 201
    ## [1] 202
    ## [1] 203
    ## [1] 204
    ## [1] 205
    ## [1] 206
    ## [1] 207
    ## [1] 208
    ## [1] 209
    ## [1] 210
    ## [1] 211
    ## [1] 212
    ## [1] 213
    ## [1] 214
    ## [1] 215
    ## [1] 216
    ## [1] 217
    ## [1] 218
    ## [1] 219
    ## [1] 220
    ## [1] 221
    ## [1] 222
    ## [1] 223
    ## [1] 224
    ## [1] 225
    ## [1] 226
    ## [1] 227
    ## [1] 228
    ## [1] 229
    ## [1] 230
    ## [1] 231
    ## [1] 232
    ## [1] 233
    ## [1] 234
    ## [1] 235
    ## [1] 236
    ## [1] 237
    ## [1] 238
    ## [1] 239
    ## [1] 240
    ## [1] 241
    ## [1] 242
    ## [1] 243
    ## [1] 244
    ## [1] 245
    ## [1] 246
    ## [1] 247
    ## [1] 248
    ## [1] 249
    ## [1] 250
    ## [1] 251
    ## [1] 252
    ## [1] 253
    ## [1] 254
    ## [1] 255
    ## [1] 256
    ## [1] 257
    ## [1] 258
    ## [1] 259
    ## [1] 260
    ## [1] 261
    ## [1] 262
    ## [1] 263
    ## [1] 264
    ## [1] 265
    ## [1] 266
    ## [1] 267
    ## [1] 268
    ## [1] 269
    ## [1] 270
    ## [1] 271
    ## [1] 272
    ## [1] 273
    ## [1] 274
    ## [1] 275
    ## [1] 276
    ## [1] 277
    ## [1] 278
    ## [1] 279
    ## [1] 280
    ## [1] 281
    ## [1] 282
    ## [1] 283
    ## [1] 284
    ## [1] 285
    ## [1] 286
    ## [1] 287
    ## [1] 288
    ## [1] 289
    ## [1] 290
    ## [1] 291
    ## [1] 292
    ## [1] 293
    ## [1] 294
    ## [1] 295
    ## [1] 296
    ## [1] 297
    ## [1] 298
    ## [1] 299
    ## [1] 300
    ## [1] 301
    ## [1] 302
    ## [1] 303
    ## [1] 304
    ## [1] 305
    ## [1] 306
    ## [1] 307
    ## [1] 308
    ## [1] 309
    ## [1] 310
    ## [1] 311
    ## [1] 312
    ## [1] 313
    ## [1] 314
    ## [1] 315
    ## [1] 316
    ## [1] 317
    ## [1] 318
    ## [1] 319
    ## [1] 320
    ## [1] 321
    ## [1] 322
    ## [1] 323
    ## [1] 324
    ## [1] 325
    ## [1] 326
    ## [1] 327
    ## [1] 328
    ## [1] 329
    ## [1] 330
    ## [1] 331
    ## [1] 332
    ## [1] 333
    ## [1] 334
    ## [1] 335
    ## [1] 336
    ## [1] 337
    ## [1] 338
    ## [1] 339
    ## [1] 340
    ## [1] 341
    ## [1] 342
    ## [1] 343
    ## [1] 344
    ## [1] 345
    ## [1] 346
    ## [1] 347
    ## [1] 348
    ## [1] 349
    ## [1] 350
    ## [1] 351
    ## [1] 352
    ## [1] 353
    ## [1] 354
    ## [1] 355
    ## [1] 356
    ## [1] 357
    ## [1] 358
    ## [1] 359
    ## [1] 360
    ## [1] 361
    ## [1] 362
    ## [1] 363
    ## [1] 364
    ## [1] 365
    ## [1] 366
    ## [1] 367
    ## [1] 368
    ## [1] 369
    ## [1] 370
    ## [1] 371
    ## [1] 372
    ## [1] 373
    ## [1] 374
    ## [1] 375
    ## [1] 376
    ## [1] 377
    ## [1] 378
    ## [1] 379
    ## [1] 380
    ## [1] 381
    ## [1] 382
    ## [1] 383
    ## [1] 384
    ## [1] 385
    ## [1] 386
    ## [1] 387
    ## [1] 388
    ## [1] 389
    ## [1] 390
    ## [1] 391
    ## [1] 392
    ## [1] 393
    ## [1] 394
    ## [1] 395
    ## [1] 396
    ## [1] 397
    ## [1] 398
    ## [1] 399
    ## [1] 400
    ## [1] 401
    ## [1] 402
    ## [1] 403
    ## [1] 404
    ## [1] 405
    ## [1] 406
    ## [1] 407
    ## [1] 408
    ## [1] 409
    ## [1] 410
    ## [1] 411
    ## [1] 412
    ## [1] 413
    ## [1] 414
    ## [1] 415
    ## [1] 416
    ## [1] 417
    ## [1] 418
    ## [1] 419
    ## [1] 420
    ## [1] 421
    ## [1] 422
    ## [1] 423
    ## [1] 424
    ## [1] 425
    ## [1] 426
    ## [1] 427
    ## [1] 428
    ## [1] 429
    ## [1] 430
    ## [1] 431
    ## [1] 432
    ## [1] 433
    ## [1] 434
    ## [1] 435
    ## [1] 436
    ## [1] 437
    ## [1] 438
    ## [1] 439
    ## [1] 440
    ## [1] 441
    ## [1] 442
    ## [1] 443
    ## [1] 444
    ## [1] 445
    ## [1] 446
    ## [1] 447
    ## [1] 448
    ## [1] 449
    ## [1] 450
    ## [1] 451
    ## [1] 452
    ## [1] 453
    ## [1] 454
    ## [1] 455
    ## [1] 456
    ## [1] 457
    ## [1] 458
    ## [1] 459
    ## [1] 460
    ## [1] 461
    ## [1] 462
    ## [1] 463
    ## [1] 464
    ## [1] 465
    ## [1] 466
    ## [1] 467
    ## [1] 468
    ## [1] 469
    ## [1] 470
    ## [1] 471
    ## [1] 472
    ## [1] 473
    ## [1] 474
    ## [1] 475
    ## [1] 476
    ## [1] 477
    ## [1] 478
    ## [1] 479
    ## [1] 480
    ## [1] 481
    ## [1] 482
    ## [1] 483
    ## [1] 484
    ## [1] 485
    ## [1] 486
    ## [1] 487
    ## [1] 488
    ## [1] 489
    ## [1] 490
    ## [1] 491
    ## [1] 492
    ## [1] 493
    ## [1] 494
    ## [1] 495
    ## [1] 496
    ## [1] 497
    ## [1] 498
    ## [1] 499
    ## [1] 500
    ## [1] 501
    ## [1] 502
    ## [1] 503
    ## [1] 504
    ## [1] 505
    ## [1] 506
    ## [1] 507
    ## [1] 508
    ## [1] 509
    ## [1] 510
    ## [1] 511
    ## [1] 512
    ## [1] 513
    ## [1] 514
    ## [1] 515
    ## [1] 516
    ## [1] 517
    ## [1] 518
    ## [1] 519
    ## [1] 520
    ## [1] 521
    ## [1] 522
    ## [1] 523
    ## [1] 524
    ## [1] 525
    ## [1] 526
    ## [1] 527
    ## [1] 528
    ## [1] 529
    ## [1] 530
    ## [1] 531
    ## [1] 532
    ## [1] 533
    ## [1] 534
    ## [1] 535
    ## [1] 536
    ## [1] 537
    ## [1] 538
    ## [1] 539
    ## [1] 540
    ## [1] 541
    ## [1] 542
    ## [1] 543
    ## [1] 544
    ## [1] 545
    ## [1] 546
    ## [1] 547
    ## [1] 548
    ## [1] 549
    ## [1] 550
    ## [1] 551
    ## [1] 552
    ## [1] 553
    ## [1] 554
    ## [1] 555
    ## [1] 556
    ## [1] 557
    ## [1] 558
    ## [1] 559
    ## [1] 560
    ## [1] 561
    ## [1] 562
    ## [1] 563
    ## [1] 564
    ## [1] 565
    ## [1] 566
    ## [1] 567
    ## [1] 568
    ## [1] 569
    ## [1] 570
    ## [1] 571
    ## [1] 572
    ## [1] 573
    ## [1] 574
    ## [1] 575
    ## [1] 576
    ## [1] 577
    ## [1] 578
    ## [1] 579
    ## [1] 580
    ## [1] 581
    ## [1] 582
    ## [1] 583
    ## [1] 584
    ## [1] 585
    ## [1] 586
    ## [1] 587
    ## [1] 588
    ## [1] 589
    ## [1] 590
    ## [1] 591
    ## [1] 592
    ## [1] 593
    ## [1] 594
    ## [1] 595
    ## [1] 596
    ## [1] 597
    ## [1] 598
    ## [1] 599
    ## [1] 600
    ## [1] 601
    ## [1] 602
    ## [1] 603
    ## [1] 604
    ## [1] 605
    ## [1] 606
    ## [1] 607
    ## [1] 608
    ## [1] 609
    ## [1] 610
    ## [1] 611
    ## [1] 612
    ## [1] 613
    ## [1] 614
    ## [1] 615
    ## [1] 616
    ## [1] 617
    ## [1] 618
    ## [1] 619
    ## [1] 620
    ## [1] 621
    ## [1] 622
    ## [1] 623
    ## [1] 624
    ## [1] 625
    ## [1] 626
    ## [1] 627
    ## [1] 628
    ## [1] 629
    ## [1] 630
    ## [1] 631
    ## [1] 632
    ## [1] 633
    ## [1] 634
    ## [1] 635
    ## [1] 636
    ## [1] 637
    ## [1] 638
    ## [1] 639
    ## [1] 640
    ## [1] 641
    ## [1] 642
    ## [1] 643
    ## [1] 644
    ## [1] 645
    ## [1] 646
    ## [1] 647
    ## [1] 648
    ## [1] 649
    ## [1] 650
    ## [1] 651
    ## [1] 652
    ## [1] 653
    ## [1] 654
    ## [1] 655
    ## [1] 656
    ## [1] 657
    ## [1] 658
    ## [1] 659
    ## [1] 660
    ## [1] 661
    ## [1] 662
    ## [1] 663
    ## [1] 664
    ## [1] 665
    ## [1] 666
    ## [1] 667
    ## [1] 668
    ## [1] 669
    ## [1] 670
    ## [1] 671
    ## [1] 672
    ## [1] 673
    ## [1] 674
    ## [1] 675
    ## [1] 676
    ## [1] 677
    ## [1] 678
    ## [1] 679
    ## [1] 680
    ## [1] 681
    ## [1] 682
    ## [1] 683
    ## [1] 684
    ## [1] 685
    ## [1] 686
    ## [1] 687
    ## [1] 688
    ## [1] 689
    ## [1] 690
    ## [1] 691
    ## [1] 692
    ## [1] 693
    ## [1] 694
    ## [1] 695
    ## [1] 696
    ## [1] 697
    ## [1] 698
    ## [1] 699
    ## [1] 700
    ## [1] 701
    ## [1] 702
    ## [1] 703
    ## [1] 704
    ## [1] 705
    ## [1] 706
    ## [1] 707
    ## [1] 708
    ## [1] 709
    ## [1] 710
    ## [1] 711
    ## [1] 712
    ## [1] 713
    ## [1] 714
    ## [1] 715
    ## [1] 716
    ## [1] 717
    ## [1] 718
    ## [1] 719
    ## [1] 720
    ## [1] 721
    ## [1] 722
    ## [1] 723
    ## [1] 724
    ## [1] 725
    ## [1] 726
    ## [1] 727
    ## [1] 728
    ## [1] 729
    ## [1] 730
    ## [1] 731
    ## [1] 732
    ## [1] 733
    ## [1] 734
    ## [1] 735
    ## [1] 736
    ## [1] 737
    ## [1] 738
    ## [1] 739
    ## [1] 740
    ## [1] 741
    ## [1] 742
    ## [1] 743
    ## [1] 744
    ## [1] 745
    ## [1] 746
    ## [1] 747
    ## [1] 748
    ## [1] 749
    ## [1] 750
    ## [1] 751
    ## [1] 752
    ## [1] 753
    ## [1] 754
    ## [1] 755
    ## [1] 756
    ## [1] 757
    ## [1] 758
    ## [1] 759
    ## [1] 760
    ## [1] 761
    ## [1] 762
    ## [1] 763
    ## [1] 764
    ## [1] 765
    ## [1] 766
    ## [1] 767
    ## [1] 768
    ## [1] 769
    ## [1] 770
    ## [1] 771
    ## [1] 772
    ## [1] 773
    ## [1] 774
    ## [1] 775
    ## [1] 776
    ## [1] 777
    ## [1] 778
    ## [1] 779
    ## [1] 780
    ## [1] 781
    ## [1] 782
    ## [1] 783
    ## [1] 784
    ## [1] 785
    ## [1] 786
    ## [1] 787
    ## [1] 788
    ## [1] 789
    ## [1] 790
    ## [1] 791
    ## [1] 792
    ## [1] 793
    ## [1] 794
    ## [1] 795
    ## [1] 796
    ## [1] 797
    ## [1] 798
    ## [1] 799
    ## [1] 800
    ## [1] 801
    ## [1] 802
    ## [1] 803
    ## [1] 804
    ## [1] 805
    ## [1] 806
    ## [1] 807
    ## [1] 808
    ## [1] 809
    ## [1] 810
    ## [1] 811
    ## [1] 812
    ## [1] 813
    ## [1] 814
    ## [1] 815
    ## [1] 816
    ## [1] 817
    ## [1] 818
    ## [1] 819
    ## [1] 820
    ## [1] 821
    ## [1] 822
    ## [1] 823
    ## [1] 824
    ## [1] 825
    ## [1] 826
    ## [1] 827
    ## [1] 828
    ## [1] 829
    ## [1] 830
    ## [1] 831
    ## [1] 832
    ## [1] 833
    ## [1] 834
    ## [1] 835
    ## [1] 836
    ## [1] 837
    ## [1] 838
    ## [1] 839
    ## [1] 840
    ## [1] 841
    ## [1] 842
    ## [1] 843
    ## [1] 844
    ## [1] 845
    ## [1] 846
    ## [1] 847
    ## [1] 848
    ## [1] 849
    ## [1] 850
    ## [1] 851
    ## [1] 852
    ## [1] 853
    ## [1] 854
    ## [1] 855
    ## [1] 856
    ## [1] 857
    ## [1] 858
    ## [1] 859
    ## [1] 860
    ## [1] 861
    ## [1] 862
    ## [1] 863
    ## [1] 864
    ## [1] 865
    ## [1] 866
    ## [1] 867
    ## [1] 868
    ## [1] 869
    ## [1] 870
    ## [1] 871
    ## [1] 872
    ## [1] 873
    ## [1] 874
    ## [1] 875
    ## [1] 876
    ## [1] 877
    ## [1] 878
    ## [1] 879
    ## [1] 880
    ## [1] 881
    ## [1] 882
    ## [1] 883
    ## [1] 884
    ## [1] 885
    ## [1] 886
    ## [1] 887
    ## [1] 888
    ## [1] 889
    ## [1] 890
    ## [1] 891
    ## [1] 892
    ## [1] 893
    ## [1] 894
    ## [1] 895
    ## [1] 896
    ## [1] 897
    ## [1] 898
    ## [1] 899
    ## [1] 900
    ## [1] 901
    ## [1] 902
    ## [1] 903
    ## [1] 904
    ## [1] 905
    ## [1] 906
    ## [1] 907
    ## [1] 908
    ## [1] 909
    ## [1] 910
    ## [1] 911
    ## [1] 912
    ## [1] 913
    ## [1] 914
    ## [1] 915
    ## [1] 916
    ## [1] 917
    ## [1] 918
    ## [1] 919
    ## [1] 920
    ## [1] 921
    ## [1] 922
    ## [1] 923
    ## [1] 924
    ## [1] 925
    ## [1] 926
    ## [1] 927
    ## [1] 928
    ## [1] 929
    ## [1] 930
    ## [1] 931
    ## [1] 932
    ## [1] 933
    ## [1] 934
    ## [1] 935
    ## [1] 936
    ## [1] 937
    ## [1] 938
    ## [1] 939
    ## [1] 940
    ## [1] 941
    ## [1] 942
    ## [1] 943
    ## [1] 944
    ## [1] 945
    ## [1] 946
    ## [1] 947
    ## [1] 948
    ## [1] 949
    ## [1] 950
    ## [1] 951
    ## [1] 952
    ## [1] 953
    ## [1] 954
    ## [1] 955
    ## [1] 956
    ## [1] 957
    ## [1] 958
    ## [1] 959
    ## [1] 960
    ## [1] 961
    ## [1] 962
    ## [1] 963
    ## [1] 964
    ## [1] 965
    ## [1] 966
    ## [1] 967
    ## [1] 968
    ## [1] 969
    ## [1] 970
    ## [1] 971
    ## [1] 972
    ## [1] 973
    ## [1] 974
    ## [1] 975
    ## [1] 976
    ## [1] 977
    ## [1] 978
    ## [1] 979
    ## [1] 980
    ## [1] 981
    ## [1] 982
    ## [1] 983
    ## [1] 984
    ## [1] 985
    ## [1] 986
    ## [1] 987
    ## [1] 988
    ## [1] 989
    ## [1] 990
    ## [1] 991
    ## [1] 992
    ## [1] 993
    ## [1] 994
    ## [1] 995
    ## [1] 996
    ## [1] 997
    ## [1] 998
    ## [1] 999
    ## [1] 1000

    head(MyObj@DyScore$Gene)

    ##   ENSMUSG00000033845   ENSMUSG00000025903   ENSMUSG00000104217 
    ##           -0.8239497           -0.8826398           -0.1752226 
    ## ENSMUSG00000025903.1 ENSMUSG00000104217.1   ENSMUSG00000033813 
    ##           -0.8826398           -0.1752226           -0.3066903

    head(MyObj@DyScore$AS)

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

    my_gtf <- system.file("extdata", "GRCm38.93_sub.gtf.gz", package = "FAScore", mustWork = TRUE)

    MyObj <- ASmapIso(MyObj, gtf = my_gtf, AStype = "exonic", cores = 10)
    MyObj <- ChooseIso(MyObj)

</br>

Adding the structural scores of transcripts

    MyObj <- matchAppris(MyObj, gtf = my_gtf, species = "MusMus")

    head(MyObj@Structure)

</br>

or other species and annotation versions downloaded from
[**APPRIS**](https://appris.bioinfo.cnio.es/#/) database

    MyObj <- matchAppris(MyObj, gtf = my_gtf, appris = OtherFile)

    head(MyObj@Structure)

</br>

Calculating the functional scores and classes

    MyObj <- CalcuFAScore(MyObj)

</br>

Sorted by FAScore

    MyObj@RFpredict <- MyObj@RFpredict[order(MyObj@RFpredict$FAScore,decreasing = T),]

    MyObj@RFpredict[,c(1,3,4,26,27)] %>% head

</br>

### Visualization

The `plotDyScore` function can rank the dynamic scores of genes and AS
events. You can label different classes of AS events, as defined by you,
along with their corresponding genes in the figure.

    ASID = list(Func = (MyObj@RFpredict$AS %>% head), nonFunc  = (MyObj@RFpredict$AS %>% tail))

    label = paste0(c(head(MyObj@RFpredict$GeneName),tail(MyObj@RFpredict$GeneName)), "|", stringr::str_split_fixed(unlist(ASID), "[|]", 2)[,2] )

    plotDyScore(MyObj, ASID = ASID, label = label, color = c( "#FF4500", "#2E9FDF"), size = 18)

</br>

The `plotFAScore` function can be used to rank the functional scores of
AS events, with the option to label specific events of interest.

    plotFAScore(MyObj, ASID = ASID, label = label, color = c( "#FF4500", "#2E9FDF"), size = 18) 

</br>  
</br>

## SessionInfo

    sessionInfo()

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
