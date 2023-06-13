Class 12 Lab: Transcriptomics and the analysis of RNA-Seq data
================
Arsam Rostami
05/12/2023

## 1. Bioconductor and DESeq2 setup

<u>**Bioconductor packages**</u> are installed differently than
“regular” R packages from CRAN. To install the core Bioconductor
packages, copy and paste the following two lines of code into your **R
console** one at a time.

``` r
# install.packages("BiocManager")
# BiocManager::install()
```

Load up the packages using the library function:

``` r
# Load BiocManager library
library(BiocManager)
```

    Bioconductor version '3.16' is out-of-date; the current release version '3.17'
      is available with R version '4.3'; see https://bioconductor.org/install

``` r
# Load DESeq2 library
library(DESeq2)
```

    Loading required package: S4Vectors

    Loading required package: stats4

    Loading required package: BiocGenerics


    Attaching package: 'BiocGenerics'

    The following objects are masked from 'package:stats':

        IQR, mad, sd, var, xtabs

    The following objects are masked from 'package:base':

        anyDuplicated, aperm, append, as.data.frame, basename, cbind,
        colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
        get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
        match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
        Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
        table, tapply, union, unique, unsplit, which.max, which.min


    Attaching package: 'S4Vectors'

    The following objects are masked from 'package:base':

        expand.grid, I, unname

    Loading required package: IRanges

    Loading required package: GenomicRanges

    Loading required package: GenomeInfoDb

    Loading required package: SummarizedExperiment

    Loading required package: MatrixGenerics

    Loading required package: matrixStats


    Attaching package: 'MatrixGenerics'

    The following objects are masked from 'package:matrixStats':

        colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
        colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
        colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
        colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
        colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
        colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
        colWeightedMeans, colWeightedMedians, colWeightedSds,
        colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
        rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
        rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
        rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
        rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
        rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
        rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
        rowWeightedSds, rowWeightedVars

    Loading required package: Biobase

    Welcome to Bioconductor

        Vignettes contain introductory material; view with
        'browseVignettes()'. To cite Bioconductor, see
        'citation("Biobase")', and for packages 'citation("pkgname")'.


    Attaching package: 'Biobase'

    The following object is masked from 'package:MatrixGenerics':

        rowMedians

    The following objects are masked from 'package:matrixStats':

        anyMissing, rowMedians

If this finished without yielding obvious error messages we can install
the **DESeq2** bioconductor package that we will use in this class:

``` r
# For this class, you'll also need DESeq2:
# BiocManager::install("DESeq2")
```

### DESeq2 Required Inputs

As input, the DESeq2 package expects (**1**) a data.frame of **count
data** (as obtained from RNA-seq or another high-throughput sequencing
experiment) and (**2**) a second data.frame with information about the
samples - often called sample metadata (or `colData` in DESeq2-speak
because it supplies metadata/information about the columns of the
countData matrix).  

# 2. Import countData and colData

Begin a new R script and use the **read.csv()** function to read these
count data and metadata files.

``` r
# Read your csv files into R
counts <- read.csv("airway_scaledcounts.csv", row.names=1)
metadata <-  read.csv("airway_metadata.csv")

# Now you can use 'counts' and 'metadata' in your analysis
```

Now, take a look at each.

``` r
# Display the first few elements counts
head(counts)
```

                    SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516
    ENSG00000000003        723        486        904        445       1170
    ENSG00000000005          0          0          0          0          0
    ENSG00000000419        467        523        616        371        582
    ENSG00000000457        347        258        364        237        318
    ENSG00000000460         96         81         73         66        118
    ENSG00000000938          0          0          1          0          2
                    SRR1039517 SRR1039520 SRR1039521
    ENSG00000000003       1097        806        604
    ENSG00000000005          0          0          0
    ENSG00000000419        781        417        509
    ENSG00000000457        447        330        324
    ENSG00000000460         94        102         74
    ENSG00000000938          0          0          0

``` r
# Display the first few rows of the data frame form metadata
head(metadata)
```

              id     dex celltype     geo_id
    1 SRR1039508 control   N61311 GSM1275862
    2 SRR1039509 treated   N61311 GSM1275863
    3 SRR1039512 control  N052611 GSM1275866
    4 SRR1039513 treated  N052611 GSM1275867
    5 SRR1039516 control  N080611 GSM1275870
    6 SRR1039517 treated  N080611 GSM1275871

Use the **View()** function to view the entire object.

``` r
# View the data frame from counts 
View(counts)
```

``` r
# View the data frame from metadata
View(metadata)
```

-   **Q1.** How many genes are in this dataset?

``` r
# Get and print the number of rows in the data frame
nrow(counts)
```

    [1] 38694

There are **38694 genes** observed using `nrow(counts)`

-   **Q2.** How many ‘control’ cell lines do we have?

``` r
# Calculate the number of 'control' in the 'dex' column
control_cell_lines<-table(metadata$dex)['control']

#print the new variable
control_cell_lines
```

    control 
          4 

There are **4 ‘control’** cell lines and 4 ‘treated’ cell lines

# 3. Toy differential gene expression

Look at the metadata object again to see which samples are `control` and
which are drug `treated`. Note that the control samples are SRR1039508,
SRR1039512, SRR1039516, and SRR1039520. This bit of code will first find
the sample `id` for those labeled control. Then calculate the mean
counts per gene across these samples:

``` r
# Get rows where 'dex' is 'control'
control <- metadata[metadata[,"dex"]=="control",]

# Get corresponding counts
control.counts <- counts[ ,control$id]

# Calculate the mean of these counts
control.mean <- rowSums( control.counts )/4 

# Display the first few elements of the calculated mean
head(control.mean)
```

    ENSG00000000003 ENSG00000000005 ENSG00000000419 ENSG00000000457 ENSG00000000460 
             900.75            0.00          520.50          339.75           97.25 
    ENSG00000000938 
               0.75 

``` r
# Filter rows where 'dex' is 'control'
metadata[,"dex"] == 'control'
```

    [1]  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE

``` r
# Filter rows where 'dex' is 'control'
metadata[metadata[,"dex"] == 'control',]
```

              id     dex celltype     geo_id
    1 SRR1039508 control   N61311 GSM1275862
    3 SRR1039512 control  N052611 GSM1275866
    5 SRR1039516 control  N080611 GSM1275870
    7 SRR1039520 control  N061011 GSM1275874

``` r
# Filter rows where 'dex' is 'control' by setting a new varibale
control <- metadata[metadata[,"dex"]=="control",]

# Access 'id' column of the 'control' data frame
control$id
```

    [1] "SRR1039508" "SRR1039512" "SRR1039516" "SRR1039520"

``` r
# Get corresponding counts
control.counts<- (counts[,control$id])

# Calculate the row means of these counts
ontrol.means <- rowMeans(control.counts)
```

An alternative way to do this same thing using the `dplyr` package from
the tidyverse is shown below.

``` r
# Install and load necessary libraries
library(dplyr)
```


    Attaching package: 'dplyr'

    The following object is masked from 'package:Biobase':

        combine

    The following object is masked from 'package:matrixStats':

        count

    The following objects are masked from 'package:GenomicRanges':

        intersect, setdiff, union

    The following object is masked from 'package:GenomeInfoDb':

        intersect

    The following objects are masked from 'package:IRanges':

        collapse, desc, intersect, setdiff, slice, union

    The following objects are masked from 'package:S4Vectors':

        first, intersect, rename, setdiff, setequal, union

    The following objects are masked from 'package:BiocGenerics':

        combine, intersect, setdiff, union

    The following objects are masked from 'package:stats':

        filter, lag

    The following objects are masked from 'package:base':

        intersect, setdiff, setequal, union

``` r
# Filter rows where 'dex' is 'control' using dplyr
control <- metadata %>% filter(dex=="control")

# Get corresponding counts using dplyr
control.counts <- counts%>% select(control$id) 

# Calculate the mean of these counts
control.mean <- rowSums(control.counts)/4

# Display the first few elements of the calculated mean
head(control.mean)
```

    ENSG00000000003 ENSG00000000005 ENSG00000000419 ENSG00000000457 ENSG00000000460 
             900.75            0.00          520.50          339.75           97.25 
    ENSG00000000938 
               0.75 

-   **Q3.** How would you make the above code in either approach more
    robust?

`rowSums(control.counts)/4` is problematic because it is not as
reproducible if new data had more or less variables.
rowMeans(control.counts) is more reproducible because it will adapt to
elements of variable added or taken away.

-   **Q4.** Follow the same procedure for the `treated` samples
    (i.e. calculate the mean per gene across drug treated samples and
    assign to a labeled vector called `treated.mean`)

``` r
# Filter rows where 'dex' is 'treated'
treated <- metadata[metadata[,"dex"]=="treated",]

# Get corresponding counts
treated.counts <- (counts[,treated$id])

# Calculate the row means of these counts
treated.mean <- rowMeans(treated.counts)

# Display the first few elements of the calculated mean
head(treated.mean)
```

    ENSG00000000003 ENSG00000000005 ENSG00000000419 ENSG00000000457 ENSG00000000460 
             658.00            0.00          546.00          316.50           78.75 
    ENSG00000000938 
               0.00 

We will combine our meancount data for bookkeeping purposes.

``` r
# Create a new data frame with 'control.mean' and 'treated.mean' as columns
meancounts <- data.frame(control.mean, treated.mean)

# Calculate the column sums of this new data frame
colSums(meancounts)
```

    control.mean treated.mean 
        23005324     22196524 

-   **Q5 (a).** Create a scatter plot showing the mean of the treated
    samples against the mean of the control samples. Your plot should
    look something like the following.

``` r
# Plot 'control.mean' against 'treated.mean'
plot(meancounts[,1],meancounts[,2], xlab="Control", ylab="Treated")
```

![](Untitled_files/figure-gfm/unnamed-chunk-19-1.png)

-   **Q5 (b).**You could also use the **ggplot2** package to make this
    figure producing the plot below. What **geom\_?()** function would
    you use for this plot?

`geom_point` function would be used for this plot

``` r
# Load necessary libraries
library(ggplot2)

# Plot 'control.mean' against 'treated.mean'
plot(meancounts[,1],meancounts[,2], xlab="Control", ylab="Treated") +
  geom_point()
```

![](Untitled_files/figure-gfm/unnamed-chunk-20-1.png)

    NULL

##### Wait a sec. There are 60,000-some rows in this data, but I’m only seeing a few dozen dots at most outside of the big clump around the origin.

-   **Q6.** Try plotting both axes on a log scale. What is the argument
    to **plot()** that allows you to do this?

Using the **`log(meancounts)`** will give desired conditions for the
plot

``` r
# Apply a log-transformation to the mean counts and create a scatter plot of the log-transformed mean counts
plot(log(meancounts))
```

![](Untitled_files/figure-gfm/unnamed-chunk-21-1.png)

There is no logarithmic expression for zeros which is omitting 15281
y-values and 15832 x-values.

### Fold Change Logarithmic

``` r
# simple division
40/20
```

    [1] 2

If there is no difference in expression then we will get zero

``` r
# logarithimic function with base 2
log2(20/20)
```

    [1] 0

If we double the expression we will get 1 fold

``` r
# logarithmic function with base 2
log2(40/20)
```

    [1] 1

these differences can be quantified with standard method with this
package `log2()`

To calculate the log2 of the fold change between the treated and
control…

``` r
# Calculate log2 fold change and add it to the 'meancounts' data frame
meancounts$log2fc <- log2(meancounts$treated.mean/
                            meancounts$control.mean)

# Display the first few elements of meancounts
head(meancounts)
```

                    control.mean treated.mean      log2fc
    ENSG00000000003       900.75       658.00 -0.45303916
    ENSG00000000005         0.00         0.00         NaN
    ENSG00000000419       520.50       546.00  0.06900279
    ENSG00000000457       339.75       316.50 -0.10226805
    ENSG00000000460        97.25        78.75 -0.30441833
    ENSG00000000938         0.75         0.00        -Inf

To remove zero values

``` r
# Identify rows with zero values in either of the first two columns
zero.vals <- which(meancounts[,1:2]==0, arr.ind=TRUE)

# Get unique row indices
to.rm <- unique(zero.vals[,1])

# Remove rows with zero values to create 'mycounts
mycounts <- meancounts[-to.rm,]

# Display the first few elements of mycounts
head(mycounts)
```

                    control.mean treated.mean      log2fc
    ENSG00000000003       900.75       658.00 -0.45303916
    ENSG00000000419       520.50       546.00  0.06900279
    ENSG00000000457       339.75       316.50 -0.10226805
    ENSG00000000460        97.25        78.75 -0.30441833
    ENSG00000000971      5219.00      6687.50  0.35769358
    ENSG00000001036      2327.00      1785.75 -0.38194109

Overexpressed and Underexpressed genes:

-   **Q7.** What is the purpose of the `arr.ind` argument in the
    **which()** function call above? Why would we then take the first
    column of the output and need to call the **unique()** function?

Let’s filter the dataset both ways to see how many genes

``` r
# Set up a new varibale for up-regulated genes higher than 2fc level pulled form mycounts data frame
up.ind <- mycounts$log2fc > 2

# Set up a new varibale for down-regulated genes lower than 2fc level pulled form mycounts data frame
down.ind <- mycounts$log2fc < (-2)
```

-   **Q8.** Using the `up.ind` vector above can you determine how many
    up regulated genes we have at the greater than 2 fc level?

``` r
# Count TRUE values in 'up.ind'
table(up.ind)['TRUE']
```

    TRUE 
     250 

Using the vector `up.ind` gave us 250 genes up-regulated

-   **Q9.** Using the `down.ind` vector above can you determine how many
    down regulated genes we have at the greater than 2 fc level?

``` r
# Count TRUE values in 'down.ind'
table(down.ind)['TRUE']
```

    TRUE 
     367 

Using the vector `down.ind` gave us 367 genes down-regulated

-   **Q10.** Do you trust these results? Why or why not?

**Yes!** As a user perspective these packages give proper statistics
using <u>negative binomial model</u> to display background expressions.
There is less fallacies made creating these rigorous data when comparing
to human problem solving if it was to be solved by hand

# 4. DESeq2 analysis

The `DESeq()` function takes a DESeqDataSet and returns a DESeqDataSet,

``` r
# Load DESeq2 library
library(DESeq2)
```

``` r
# Use citation function to pull relavent data regarding DESeq2 within published work
citation("DESeq2")
```


    To cite package 'DESeq2' in publications use:

      Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change
      and dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550
      (2014)

    A BibTeX entry for LaTeX users is

      @Article{,
        title = {Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2},
        author = {Michael I. Love and Wolfgang Huber and Simon Anders},
        year = {2014},
        journal = {Genome Biology},
        doi = {10.1186/s13059-014-0550-8},
        volume = {15},
        issue = {12},
        pages = {550},
      }

#### Create specific objects to apply the functions creating an object for the whole experiment with design in mind (treated vs control)

Let’s generate the specific object that DESeq2 needs:

``` r
# Set a new variable as dds to build the required DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData=counts, 
                              colData=metadata, 
                              design=~dex)
```

    converting counts to integer mode

    Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    design formula are characters, converting to factors

``` r
# Print dds
dds
```

    class: DESeqDataSet 
    dim: 38694 8 
    metadata(1): version
    assays(1): counts
    rownames(38694): ENSG00000000003 ENSG00000000005 ... ENSG00000283120
      ENSG00000283123
    rowData names(0):
    colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
    colData names(4): id dex celltype geo_id

## DESeq analysis

Would result an error “couldn’t find results. you should first run
DESeq()”

``` r
# results(dds)
```

DESeq pipeline on the `dds` object, and reassigning the whole thing back
to `dds`

``` r
# set the DESeq() fucntion using dds to a new dds variable 
dds <- DESeq(dds)
```

    estimating size factors

    estimating dispersions

    gene-wise dispersion estimates

    mean-dispersion relationship

    final dispersion estimates

    fitting model and testing

``` r
# print dds
dds
```

    class: DESeqDataSet 
    dim: 38694 8 
    metadata(1): version
    assays(4): counts mu H cooks
    rownames(38694): ENSG00000000003 ENSG00000000005 ... ENSG00000283120
      ENSG00000283123
    rowData names(22): baseMean baseVar ... deviance maxCooks
    colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
    colData names(5): id dex celltype geo_id sizeFactor

### Getting results

we can get results out of the object simply by calling
the `results()` function

``` r
# Run the DESeq pipeline adjusting the previous data
res <- results(dds)

# prit res
res
```

    log2 fold change (MLE): dex treated vs control 
    Wald test p-value: dex treated vs control 
    DataFrame with 38694 rows and 6 columns
                     baseMean log2FoldChange     lfcSE      stat    pvalue
                    <numeric>      <numeric> <numeric> <numeric> <numeric>
    ENSG00000000003  747.1942     -0.3507030  0.168246 -2.084470 0.0371175
    ENSG00000000005    0.0000             NA        NA        NA        NA
    ENSG00000000419  520.1342      0.2061078  0.101059  2.039475 0.0414026
    ENSG00000000457  322.6648      0.0245269  0.145145  0.168982 0.8658106
    ENSG00000000460   87.6826     -0.1471420  0.257007 -0.572521 0.5669691
    ...                   ...            ...       ...       ...       ...
    ENSG00000283115  0.000000             NA        NA        NA        NA
    ENSG00000283116  0.000000             NA        NA        NA        NA
    ENSG00000283119  0.000000             NA        NA        NA        NA
    ENSG00000283120  0.974916      -0.668258   1.69456 -0.394354  0.693319
    ENSG00000283123  0.000000             NA        NA        NA        NA
                         padj
                    <numeric>
    ENSG00000000003  0.163035
    ENSG00000000005        NA
    ENSG00000000419  0.176032
    ENSG00000000457  0.961694
    ENSG00000000460  0.815849
    ...                   ...
    ENSG00000283115        NA
    ENSG00000283116        NA
    ENSG00000283119        NA
    ENSG00000283120        NA
    ENSG00000283123        NA

We can summarize some basic tallies using the summary function.

``` r
# Use summary function to achieve p-value at 0.05 for the argument 
summary(res, alpha=0.05)
```


    out of 25258 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 1242, 4.9%
    LFC < 0 (down)     : 939, 3.7%
    outliers [1]       : 142, 0.56%
    low counts [2]     : 9971, 39%
    (mean count < 10)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

  
The results function contains a number of arguments to customize the
results table. By default the argument `alpha` is set to 0.1.

``` r
# Extract results using p-value of 0.05 as the parameter
res05 <- results(dds, alpha=0.05)

# Print the summary
summary(res05)
```


    out of 25258 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 1236, 4.9%
    LFC < 0 (down)     : 933, 3.7%
    outliers [1]       : 142, 0.56%
    low counts [2]     : 9033, 36%
    (mean count < 6)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

# 5. Adding annotation data

Here we load the **AnnotationDbi** package and the annotation data
package for humans **org.Hs.eg.db**.

``` r
# Load AnnotationDbi library
library("AnnotationDbi")
```


    Attaching package: 'AnnotationDbi'

    The following object is masked from 'package:dplyr':

        select

``` r
# Load org.Hs.eg.db library
library("org.Hs.eg.db")
```

To get a list of all available key types that we can use to map between,
use the `columns()` function:

``` r
# Get the list of available columns or key types in the org.Hs.eg.db packag
columns(org.Hs.eg.db)
```

     [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
     [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
    [11] "GENETYPE"     "GO"           "GOALL"        "IPI"          "MAP"         
    [16] "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"        
    [21] "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
    [26] "UNIPROT"     

The main function we will use from the **AnnotationDbi** package is
called **mapIds()**. We provide the row names of our results table as a
key, and specify that `keytype=ENSEMBL.`

``` r
# Add gene symbols to the results
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res), # Our genenames
                     keytype="ENSEMBL",        # The format of our genenames
                     column="SYMBOL",          # The new format we want to add
                     multiVals="first")
```

    'select()' returned 1:many mapping between keys and columns

``` r
# Display the first few elements of res
head(res)
```

    log2 fold change (MLE): dex treated vs control 
    Wald test p-value: dex treated vs control 
    DataFrame with 6 rows and 7 columns
                      baseMean log2FoldChange     lfcSE      stat    pvalue
                     <numeric>      <numeric> <numeric> <numeric> <numeric>
    ENSG00000000003 747.194195     -0.3507030  0.168246 -2.084470 0.0371175
    ENSG00000000005   0.000000             NA        NA        NA        NA
    ENSG00000000419 520.134160      0.2061078  0.101059  2.039475 0.0414026
    ENSG00000000457 322.664844      0.0245269  0.145145  0.168982 0.8658106
    ENSG00000000460  87.682625     -0.1471420  0.257007 -0.572521 0.5669691
    ENSG00000000938   0.319167     -1.7322890  3.493601 -0.495846 0.6200029
                         padj      symbol
                    <numeric> <character>
    ENSG00000000003  0.163035      TSPAN6
    ENSG00000000005        NA        TNMD
    ENSG00000000419  0.176032        DPM1
    ENSG00000000457  0.961694       SCYL3
    ENSG00000000460  0.815849    C1orf112
    ENSG00000000938        NA         FGR

-   **Q11.** Run the **mapIds()** function two more times to add the
    Entrez ID and UniProt accession and GENENAME as new columns called
    `res$entrez`, `res$uniprot` and `res$genename`.

``` r
# Add Entrez ID to the results
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
```

    'select()' returned 1:many mapping between keys and columns

``` r
# Add UniProt accession to the results
res$uniprot <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="UNIPROT",
                     keytype="ENSEMBL",
                     multiVals="first")
```

    'select()' returned 1:many mapping between keys and columns

``` r
# Add gene names to the results
res$genename <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="GENENAME",
                     keytype="ENSEMBL",
                     multiVals="first")
```

    'select()' returned 1:many mapping between keys and columns

``` r
# Note: GENENAME column is not directly available in org.Hs.eg.db.
# Therefore, use SYMBOL column as an alternative.

# Print the updated 'res' data frame
head(res)
```

    log2 fold change (MLE): dex treated vs control 
    Wald test p-value: dex treated vs control 
    DataFrame with 6 rows and 10 columns
                      baseMean log2FoldChange     lfcSE      stat    pvalue
                     <numeric>      <numeric> <numeric> <numeric> <numeric>
    ENSG00000000003 747.194195     -0.3507030  0.168246 -2.084470 0.0371175
    ENSG00000000005   0.000000             NA        NA        NA        NA
    ENSG00000000419 520.134160      0.2061078  0.101059  2.039475 0.0414026
    ENSG00000000457 322.664844      0.0245269  0.145145  0.168982 0.8658106
    ENSG00000000460  87.682625     -0.1471420  0.257007 -0.572521 0.5669691
    ENSG00000000938   0.319167     -1.7322890  3.493601 -0.495846 0.6200029
                         padj      symbol      entrez     uniprot
                    <numeric> <character> <character> <character>
    ENSG00000000003  0.163035      TSPAN6        7105  A0A024RCI0
    ENSG00000000005        NA        TNMD       64102      Q9H2S6
    ENSG00000000419  0.176032        DPM1        8813      O60762
    ENSG00000000457  0.961694       SCYL3       57147      Q8IZE3
    ENSG00000000460  0.815849    C1orf112       55732  A0A024R922
    ENSG00000000938        NA         FGR        2268      P09769
                                  genename
                               <character>
    ENSG00000000003          tetraspanin 6
    ENSG00000000005            tenomodulin
    ENSG00000000419 dolichyl-phosphate m..
    ENSG00000000457 SCY1 like pseudokina..
    ENSG00000000460 chromosome 1 open re..
    ENSG00000000938 FGR proto-oncogene, ..

You can arrange and view the results by the adjusted p-value

``` r
# set a new varibale using order function 
ord <- order( res$padj )

#View(res[ord,])
head(res[ord,])
```

    log2 fold change (MLE): dex treated vs control 
    Wald test p-value: dex treated vs control 
    DataFrame with 6 rows and 10 columns
                     baseMean log2FoldChange     lfcSE      stat      pvalue
                    <numeric>      <numeric> <numeric> <numeric>   <numeric>
    ENSG00000152583   954.771        4.36836 0.2371268   18.4220 8.74490e-76
    ENSG00000179094   743.253        2.86389 0.1755693   16.3120 8.10784e-60
    ENSG00000116584  2277.913       -1.03470 0.0650984  -15.8944 6.92855e-57
    ENSG00000189221  2383.754        3.34154 0.2124058   15.7319 9.14433e-56
    ENSG00000120129  3440.704        2.96521 0.2036951   14.5571 5.26424e-48
    ENSG00000148175 13493.920        1.42717 0.1003890   14.2164 7.25128e-46
                           padj      symbol      entrez     uniprot
                      <numeric> <character> <character> <character>
    ENSG00000152583 1.32441e-71     SPARCL1        8404  A0A024RDE1
    ENSG00000179094 6.13966e-56        PER1        5187      O15534
    ENSG00000116584 3.49776e-53     ARHGEF2        9181      Q92974
    ENSG00000189221 3.46227e-52        MAOA        4128      P21397
    ENSG00000120129 1.59454e-44       DUSP1        1843      B4DU40
    ENSG00000148175 1.83034e-42        STOM        2040      F8VSL7
                                  genename
                               <character>
    ENSG00000152583           SPARC like 1
    ENSG00000179094 period circadian reg..
    ENSG00000116584 Rho/Rac guanine nucl..
    ENSG00000189221    monoamine oxidase A
    ENSG00000120129 dual specificity pho..
    ENSG00000148175               stomatin

Finally, let’s write out the ordered significant results with
annotations.

``` r
# Save the ordered results to a CSV file
write.csv(res[ord,], "deseq_results.csv")
```

# 6. Data Visualization

## Volcano plots

``` r
# Assuming 'res' is the results object obtained from DESeq analysis

# Create the scatter plot
plot( res$log2FoldChange,  -log10(res$padj), 
      xlab="Log2(FoldChange)",
      ylab="-Log10(P-value)")
```

![](Untitled_files/figure-gfm/unnamed-chunk-45-1.png)

To make this more useful we can add some guidelines (with the `abline()`
function)

``` r
# Create the scatter plot
plot( res$log2FoldChange,  -log10(res$padj), 
 ylab="-Log(P-value)", xlab="Log2(FoldChange)")

# Add some cut-off lines
abline(v=c(-2,2), col="darkgray", lty=2)
abline(h=-log(0.05), col="darkgray", lty=2)
```

![](Untitled_files/figure-gfm/unnamed-chunk-46-1.png)

To color the points we will setup a custom color vector indicating
transcripts with large fold change and significant differences between
conditions:

``` r
# Setup our custom point color vector 
mycols <- rep("gray", nrow(res))
mycols[ abs(res$log2FoldChange) > 2 ]  <- "red" 

inds <- (res$padj < 0.01) & (abs(res$log2FoldChange) > 2 )
mycols[ inds ] <- "blue"

# Volcano plot with custom colors 
plot( res$log2FoldChange,  -log10(res$padj), 
 col=mycols, ylab="-Log(P-value)", xlab="Log2(FoldChange)" )

# Cut-off lines
abline(v=c(-2,2), col="gray", lty=2)
abline(h=-log(0.1), col="gray", lty=2)
```

![](Untitled_files/figure-gfm/unnamed-chunk-47-1.png)

``` r
# BiocManager::install("EnhancedVolcano")
```

``` r
# Load up neccessary libraries 
library(EnhancedVolcano)
```

    Loading required package: ggrepel

``` r
# Convert 'res' to a data frame
x<- as.data.frame(res)
```

``` r
# Create the enhanced volcano plot
EnhancedVolcano(x,
    lab = x$symbol,
    x = 'log2FoldChange',
    y = 'pvalue')
```

![](Untitled_files/figure-gfm/unnamed-chunk-51-1.png)

# 7. Pathway analysis

First we need to do our one time install of these required bioconductor
packages:

``` r
# Run in your R console (i.e. not your Rmarkdown doc!)
# BiocManager::install( c("pathview", "gage", "gageData") )
```

``` r
# Load necessary libraries
library(pathview)
```

    ##############################################################################
    Pathview is an open source software package distributed under GNU General
    Public License version 3 (GPLv3). Details of GPLv3 is available at
    http://www.gnu.org/licenses/gpl-3.0.html. Particullary, users are required to
    formally cite the original Pathview paper (not just mention it) in publications
    or products. For details, do citation("pathview") within R.

    The pathview downloads and uses KEGG data. Non-academic uses may require a KEGG
    license agreement (details at http://www.kegg.jp/kegg/legal.html).
    ##############################################################################

``` r
library(gage)
```

``` r
library(gageData)

# Load KEGG pathway sets for humans
data(kegg.sets.hs)

# Examine the first 2 pathways in this kegg set for humans
head(kegg.sets.hs, 2)
```

    $`hsa00232 Caffeine metabolism`
    [1] "10"   "1544" "1548" "1549" "1553" "7498" "9"   

    $`hsa00983 Drug metabolism - other enzymes`
     [1] "10"     "1066"   "10720"  "10941"  "151531" "1548"   "1549"   "1551"  
     [9] "1553"   "1576"   "1577"   "1806"   "1807"   "1890"   "221223" "2990"  
    [17] "3251"   "3614"   "3615"   "3704"   "51733"  "54490"  "54575"  "54576" 
    [25] "54577"  "54578"  "54579"  "54600"  "54657"  "54658"  "54659"  "54963" 
    [33] "574537" "64816"  "7083"   "7084"   "7172"   "7363"   "7364"   "7365"  
    [41] "7366"   "7367"   "7371"   "7372"   "7378"   "7498"   "79799"  "83549" 
    [49] "8824"   "8833"   "9"      "978"   

Note that we used the **mapIDs()** function above to obtain Entrez gene
IDs (stored in `res$entrez`) and we have the fold change results from
DESeq2 analysis (stored in `res$log2FoldChange`).

``` r
# Extract log2 fold changes from 'res'
foldchanges = res$log2FoldChange

# Assign Entrez gene IDs as names to the fold changes
names(foldchanges) = res$entrez

# Display the first few entries of fold changes
head(foldchanges)
```

           7105       64102        8813       57147       55732        2268 
    -0.35070302          NA  0.20610777  0.02452695 -0.14714205 -1.73228897 

Now, let’s run the **gage** pathway analysis.

``` r
# Get the results
keggres = gage(foldchanges, gsets=kegg.sets.hs)
```

Now lets look at the object returned from **gage()**.

``` r
# View the attributes of 'keggres'
attributes(keggres)
```

    $names
    [1] "greater" "less"    "stats"  

Like any list we can use the dollar syntax to access a named element,
e.g. `head(keggres$greater)` and`head(keggres$less)`.

Lets look at the first few down (less) pathway results:

``` r
# Look at the first three down (less) pathways
head(keggres$less, 3)
```

                                          p.geomean stat.mean        p.val
    hsa05332 Graft-versus-host disease 0.0004250461 -3.473346 0.0004250461
    hsa04940 Type I diabetes mellitus  0.0017820293 -3.002352 0.0017820293
    hsa05310 Asthma                    0.0020045888 -3.009050 0.0020045888
                                            q.val set.size         exp1
    hsa05332 Graft-versus-host disease 0.09053483       40 0.0004250461
    hsa04940 Type I diabetes mellitus  0.14232581       42 0.0017820293
    hsa05310 Asthma                    0.14232581       29 0.0020045888

To begin with lets manually supply a `pathway.id` (namely the first part
of the `"hsa05310 Asthma"`) that we could see from the print out above.

``` r
# Load necessary libraries
library(pathview)

# Generate pathway visualization
pathview(gene.data=foldchanges, pathway.id="hsa05310")
```

    'select()' returned 1:1 mapping between keys and columns

    Info: Working in directory /Users/arsamrostami/Desktop/bimm143/Class 12 Lab

    Info: Writing image file hsa05310.pathview.png

![![](hsa05310.pathview.png)](hsa05310.png)

``` r
# A different PDF based output of the same data
pathview(gene.data=foldchanges, pathway.id="hsa05310", kegg.native=FALSE)
```

    'select()' returned 1:1 mapping between keys and columns

    Info: Working in directory /Users/arsamrostami/Desktop/bimm143/Class 12 Lab

    Info: Writing image file hsa05310.pathview.pdf

**Q12**. Can you do the same procedure as above to plot the pathview
figures for the top 2 down-reguled pathways?
