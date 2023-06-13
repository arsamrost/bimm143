Class 18 Mini-Project: Mutational Signatures in Human Cancer
================
Arsam Rosatmi
06/01/2023

# Background

it is possible to **identify the causes leading to a specific tumor.**
They are the consequence of multiple mutational processes, which
generate unique combinations of mutation types, termed **mutational
signatures**.

# 1. Exploring a cancer sequencing data portal

Our input data will be cancer sequencing data from the publicly
available [cBioPortal](https://www.cbioportal.org/) platform

Group 1 will analyze **liver hepatocellular carcinomas**, group 2 **lung
squamous cell carcinomas**, and group 3 **skin cutaneous melanoma**.

First of all, after selecting the data set of choice, we will choose
the `Explore selected studies` button. This will lead us to
the `Summary panel`, where different clinical features can be explored.

### **Discussion \#1 (Lung Cancer Data)**

Please work with your group to explore the clinical and molecular data
categories present in cBioPortal and answer the questions below.

**Q.** How many cancer samples are included in the dataset?

-   487 samples for **Lung Squamous Cell Carcinoma**

**Q.** Which is the most mutated gene?

-   TP53

**Q.** Which is the most common treatment undergone by patients?

-   Cisplatin,

# 2. Downloading cancer sequencing data- Lung

To download the cancer sequencing data to use for the mutational
signature analysis, we will use the following button within
the `Summary` panel. This will download all available data for the
cancer data set, including clinical and molecular data. The file
containing the mutation data is `data_mutations.txt`. 

# 3. Generating mutational matrices and visualizing mutational profiles

Benefiting from the standard MAF format, we can use the `maftools` R
package to manage this sequencing data. **Following the code below, you
can install the `maftools` R package, as well as read your input data
and generate mutational matrices.** 

In this case, we are going to focus on the **SBS96 mutational context**,
which, as mentioned in the lecture, allows classifying single nucleotide
mutations in different categories based on the mutated nucleotide

``` r
# Install required packages
if (!require("BiocManager")){
    install.packages("BiocManager")
}
```

    Loading required package: BiocManager

    Bioconductor version '3.16' is out-of-date; the current release version '3.17'
      is available with R version '4.3'; see https://bioconductor.org/install

### More Packages…

``` r
# Install required packages
if (!require("maftools")){
    BiocManager::install("maftools")
}
```

    Loading required package: maftools

### And more…

``` r
# Install required packages
if (!require("BSgenome.Hsapiens.UCSC.hg19")){         # reference genome needed to
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")   # generate mutational matrices
}
```

    Loading required package: BSgenome.Hsapiens.UCSC.hg19

    Loading required package: BSgenome

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

    Loading required package: S4Vectors

    Loading required package: stats4


    Attaching package: 'S4Vectors'

    The following objects are masked from 'package:base':

        expand.grid, I, unname

    Loading required package: IRanges

    Loading required package: GenomeInfoDb

    Loading required package: GenomicRanges

    Loading required package: Biostrings

    Loading required package: XVector


    Attaching package: 'Biostrings'

    The following object is masked from 'package:base':

        strsplit

    Loading required package: rtracklayer

For the **visualization of SBS96 mutational profiles**, we will make use
of the `MutationalPatterns` R package. This library is commonly used for
all kinds of mutational signature analysis, and we will also use it for
the subsequent assignment analysis.

``` r
# Read maf file
library(maftools)
coad = read.maf('lusc_tcga_pan_can_atlas_2018/data_mutations.txt')
```

    -Reading
    -Validating
    --Removed 14038 duplicated variants
    -Silent variants: 60539 
    -Summarizing
    --Possible FLAGS among top ten genes:
      TTN
      MUC16
      USH2A
      SYNE1
    -Processing clinical data
    --Missing clinical data
    -Finished in 7.613s elapsed (6.893s cpu) 

``` r
# Generate mutational matrix (SBS96 context)
mm_coad = trinucleotideMatrix(maf = coad, prefix = 'chr', add = TRUE,
                              ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
```

    -Extracting 5' and 3' adjacent bases
    -Extracting +/- 20bp around mutated bases for background C>T estimation
    -Estimating APOBEC enrichment scores
    --Performing one-way Fisher's test for APOBEC enrichment
    ---APOBEC related mutations are enriched in  30.128 % of samples (APOBEC enrichment score > 2 ;  141  of  468  samples)
    -Creating mutation matrix
    --matrix of dimension 469x96

``` r
mm_coad = t(mm_coad$nmf_matrix)
```

### Additional packages:

``` r
# Install MutationalPatterns package
if (!require("MutationalPatterns")){
BiocManager::install("MutationalPatterns")
}
```

    Loading required package: MutationalPatterns

    Loading required package: NMF

    Loading required package: registry

    Loading required package: rngtools

    Loading required package: cluster

    NMF - BioConductor layer [OK] | Shared memory capabilities [NO: bigmemory] | Cores 2/2

      To enable shared memory capabilities, try: install.extras('
    NMF
    ')


    Attaching package: 'NMF'

    The following object is masked from 'package:S4Vectors':

        nrun

``` r
# Generate mutational profiles (4 random samples)
library(MutationalPatterns)
set.seed(11111) # fixing the seed for random number generation

# Generate plot for each nucleoide change
samples_to_plot = sample(1:ncol(mm_coad),4) # selecting 4 random samples

# Plot the 96 profile
plot_96_profile(mm_coad[,samples_to_plot], condensed = T)
```

![](Class-18-Cancer_files/figure-gfm/unnamed-chunk-6-1.png)

Plots the samples for mutations:

``` r
# Generate mutational profiles (top 4 mutated samples and top 4 less mutated)
# using colsums() to generate all the columns in summation
mutations_in_samples = colSums(mm_coad)

# use the sort function in a decreasing order using sort()
mutations_in_samples = sort(mutations_in_samples, decreasing = T)

# Generate a plot of the top four neucleotide changes
samples_to_plot = names(mutations_in_samples)[1:4]

# Plot the 96 profile
plot_96_profile(mm_coad[,samples_to_plot], condensed = T)
```

![](Class-18-Cancer_files/figure-gfm/unnamed-chunk-7-1.png)

``` r
# Sort the mutations
mutations_in_samples = sort(mutations_in_samples, decreasing = F)

# Get the samples to plot
samples_to_plot = names(mutations_in_samples)[1:4]

# Plot the 96 profile
plot_96_profile(mm_coad[,samples_to_plot], condensed = T)
```

![](Class-18-Cancer_files/figure-gfm/unnamed-chunk-8-1.png)

``` r
# Generate average mutational profiles
relative_mutational_profile = apply(mm_coad, 2, prop.table) # obtained relative
                                                            #mutational matrix

# Calculate the average mutational profile
average_mutational_profile = rowMeans(relative_mutational_profile)

# Convert the average profile to a data frame
average_mutational_profile = data.frame(average_mutational_profile)

# Plot the 96 profile
plot_96_profile(average_mutational_profile, condensed = T)
```

![](Class-18-Cancer_files/figure-gfm/unnamed-chunk-9-1.png)

# 4. COSMIC reference mutational signatures

These reference mutational signatures have been associated with
specific **environmental exposures, lifestyle choices, and endogenous
cellular mechanisms**. 

![](images/Screenshot%202023-06-02%20at%202.38.04%20PM.png)

# 5. Assigning reference mutational signatures

First, a standard assignment analysis carried out by
the`fit_to_signatures` function using the non-negative least squares
(NNLS) algorithm. On the other hand, we will also use
the`fit_to_signatures_strict` function.

``` r
# Mutational signature assignment
cosmic_signatures = get_known_signatures(source = 'COSMIC_v3.2')
fit_res = fit_to_signatures(mm_coad, cosmic_signatures)

# Top contributing signatures
contributions = fit_res$contribution

# Calculate the means of the rows
top_contributing_signatures_abs = rowMeans(contributions)

# Sort the values in decreasing order and keep only the top 4
top_contributing_signatures_abs = sort(top_contributing_signatures_abs,
                                       decreasing = T)[1:4]

# Top 4 contributing signatures (absolute values)
top_contributing_signatures_abs
```

        SBS4    SBS24    SBS39    SBS13 
    104.9223  30.2159  29.8143  26.5996 

-   SBS4: Associated with tobacco smoking. Its profile is similar to the
    mutational spectrum observed in experimental systems exposed to
    tobacco carcinogens such as benzo\[a\]pyrene. SBS4 is, therefore,
    likely due to direct DNA damage by tobacco smoke mutagens.

-   SBS24: Aflatoxin exposure. SBS24 has been found in cancer samples
    with known exposures to aflatoxin and the pattern of mutations
    exhibited by the signature is consistent with that observed in
    experimental systems exposed to aflatoxin.

-   SBS39: Unknown.

-   SBS13: Attributed to activity of the AID/APOBEC family of cytidine
    deaminases on the basis of similarities in the sequence context of
    cytosine mutations caused by APOBEC enzymes in experimental systems.
    APOBEC3A is probably responsible for most mutations in human cancer,
    although APOBEC3B may also contribute (these differ in the sequence
    context two bases 5' to the mutated cytosine, see 1536 mutation
    classification signature extraction). SBS13 mutations are likely
    generated by error prone polymerases (such as REV1) replicating
    across abasic sites generated by base excision repair removal of
    uracil.

to get the relative contribution for the mutation signatures:

``` r
# Calculate relative contributions
relative_contributions = apply(contributions,2,prop.table)

# Calculate the means of the rows
top_contributing_signatures_rel = rowMeans(relative_contributions)

# Sort the values in decreasing order and keep only the top 4
top_contributing_signatures_rel = sort(top_contributing_signatures_rel,
                                       decreasing = T)[1:4]

## Top 4 contributing signatures (relative values)
top_contributing_signatures_rel
```

          SBS4      SBS24      SBS39      SBS13 
    0.23615327 0.08534249 0.08064435 0.06469107 

``` r
# Mutational signature assignment strict
fit_res_strict = fit_to_signatures_strict(mm_coad, cosmic_signatures)

# Extract the fit results
fit_res_strict = fit_res_strict$fit_res

# Get the contributions from the fit results
contributions_strict = fit_res_strict$contribution
```

# 6. Visualizing mutational signature assignment results

``` r
# Visualization of signature assignment results (fit_to_signatures)
set.seed(11111)

# Sample 4 columns from the data
samples_to_plot = sample(1:ncol(mm_coad),4)

# Plot the contributions for the selected samples
plot_contribution(contributions[,samples_to_plot], mode = "absolute")
```

![](Class-18-Cancer_files/figure-gfm/unnamed-chunk-13-1.png)

``` r
# Plot the contributions for the selected samples
plot_contribution(contributions[,samples_to_plot], mode = "relative")
```

![](Class-18-Cancer_files/figure-gfm/unnamed-chunk-14-1.png)

``` r
# Plot the contributions using a heatmap, without clustering the samples
plot_contribution_heatmap(contributions, cluster_samples = F)
```

![](Class-18-Cancer_files/figure-gfm/unnamed-chunk-15-1.png)

``` r
# Visualization of signature assignment results (strict)
plot_contribution(contributions_strict[,samples_to_plot], mode = "absolute")
```

![](Class-18-Cancer_files/figure-gfm/unnamed-chunk-16-1.png)

``` r
# Visualization of signature assignment results (strict)
plot_contribution(contributions_strict[,samples_to_plot], mode = "relative")
```

![](Class-18-Cancer_files/figure-gfm/unnamed-chunk-17-1.png)

``` r
# Plot the contributions using a heatmap, without clustering the samples
plot_contribution_heatmap(contributions_strict, cluster_samples = F)
```

![](Class-18-Cancer_files/figure-gfm/unnamed-chunk-18-1.png)

``` r
# Cosine similarity reconstruction vs. original mutational profile (fit_to_signatures)
set.seed(11111)

# Sample 4 columns from the data
samples_to_plot = sample(1:ncol(mm_coad),4)

# Plot the original vs reconstructed data for the selected samples
plot_original_vs_reconstructed(mm_coad[,samples_to_plot],
                               fit_res$reconstructed[,samples_to_plot], 
                               y_intercept = 0.90)
```

![](Class-18-Cancer_files/figure-gfm/unnamed-chunk-19-1.png)

``` r
# Cosine similarity reconstruction vs. original mutational profile (strict)
plot_original_vs_reconstructed(mm_coad[,samples_to_plot],
                               fit_res_strict$reconstructed[,samples_to_plot], 
                               y_intercept = 0.90)
```

![](Class-18-Cancer_files/figure-gfm/unnamed-chunk-20-1.png)

### **Discussion \#2**

**Q.** Which is the etiology of the top absolute contributing signature
for liver cancer?

-   Aristolochic Acid Exposure

**Q.** Which is the most prominent mutational context for the top
contributing signature in skin cancer?  

-   C\>T 

**Q.** The etiology of the top contributing signature for lung cancer
corresponds to an endogenous cellular mechanism. 

-   FALSE

**Q.** SBS4 is one of the most common signatures found in lung cancer
and is associated with tobacco smoking. 

-   TRUE

**Q.** SBS7d is one of the most common signatures in skin cancer and is
associated with UV light exposure and high numbers of C\>T mutations.

-   FALSE
