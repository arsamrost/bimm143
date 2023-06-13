Class 9: Structural Bioinformatics (Pt. 1)
================
Arsam Rostami

## PDB statistics

to read the file we are going to use the `read.csv` command and
importing the downloaded file named “Data Export Summary.csv”. Download
a CSV file from the PDB site (accessible from **“Analyze” \> “PDB
Statistics” \> “by Experimental Method and Molecular Type”**

``` r
# Read in the CSV file and assign it to the variable 'pdb_stats'
pdb_stats<-read.csv("Data Export Summary.csv", row.names=1)

# Print the first few rows of the data to verify that it was loaded correctly
pdb_stats
```

                              X.ray     EM    NMR Multiple.methods Neutron Other
    Protein (only)          154,766 10,155 12,187              191      72    32
    Protein/Oligosaccharide   9,083  1,802     32                7       1     0
    Protein/NA                8,110  3,176    283                6       0     0
    Nucleic acid (only)       2,664     94  1,450               12       2     1
    Other                       163      9     32                0       0     0
    Oligosaccharide (only)       11      0      6                1       0     4
                              Total
    Protein (only)          177,403
    Protein/Oligosaccharide  10,925
    Protein/NA               11,575
    Nucleic acid (only)       4,223
    Other                       204
    Oligosaccharide (only)       22

-   **Note**: there is 6 observation out of 7 variables

**Q1: What percentage of structures in the PDB are solved by X-Ray and
Electron Microscopy.**

I need to sum all the elements of the X.ray column.

``` r
#  Pull the X.ray column and convert numbers with commas expressed as characters to numeric values
pdb_stats$X.ray
```

    [1] "154,766" "9,083"   "8,110"   "2,664"   "163"     "11"     

``` r
# use gsub() command by having a pattern (',') and a replacement of ('') which indicates to removing the comma enabling characters ro be recognized as numerics...

gsub(',', '', pdb_stats$X.ray)
```

    [1] "154766" "9083"   "8110"   "2664"   "163"    "11"    

We are going to use `as.numeric()` command to define the numeric as a
vector

``` r
#convert the vector to numeric without paranteses by using as.numeric() command 

as.numeric(sub(',', '', pdb_stats$X.ray))
```

    [1] 154766   9083   8110   2664    163     11

I use the sum command now after removing all the varibales

``` r
# Calculate the sum of the X.ray column, converting numbers with commas expressed as characters to numeric values
n_xrays<- sum(as.numeric(sub(',', '', pdb_stats$X.ray)))

# Print the result
n_xrays
```

    [1] 174797

``` r
# Calculate the sum of the EM column, converting numbers with commas expressed as characters to numeric values
n_em<- sum(as.numeric(sub(',', '', pdb_stats$EM)))

# Print the result
n_em
```

    [1] 15236

``` r
# Calculate the sum of the Total column, converting numbers with commas expressed as characters to numeric values
n_total<- sum(as.numeric(sub(',', '', pdb_stats$Total)))

# Print the result
n_total
```

    [1] 204352

``` r
# Calculate the proportions of X-rays and EM values in the pdb_stats data
p_xrays<-(n_xrays)/ n_total
p_em<-(n_em)/ n_total


# Calculate the total proportion
p_total= p_xrays+p_em

# Print the results
p_xrays
```

    [1] 0.8553721

``` r
p_em
```

    [1] 0.07455763

``` r
p_total
```

    [1] 0.9299297

-   **92.9% s**tructures in the PDB are solved by X-Ray and Electron
    Microscopy

-   **85.5%** **s**tructures in the PDB are solved by X-Ray

-   **7.4% s**tructures in the PDB are solved by Electron Microscopy

**Q2: What proportion of structures in the PDB are protein?**

we are going to use the first row for all proteins

``` r
# using the gsub command 
# Calculate the proportion of total protein in the first row of pdb_stats, where numbers with commas expressed as characters are converted to numeric values
total_protein<-as.numeric(gsub(',', '', pdb_stats[1,7]))
```

``` r
# use divion command to achieve proportion of structures in the PDB are protein
total_protein/n_total
```

    [1] 0.8681246

-   The proportion of structures in the PDB are protein is 0.868 which
    in percentage translates to 86.8%

**Q3: Type HIV in the PDB website search box on the home page and
determine how many HIV-1 protease structures are in the current PDB?**

-   It was difficult to determine how many HIV-1 protease structures are
    in the current PDB

## The PDB format

# 2. Visualizing the HIV-1 protease structure

**Q4: Water molecules normally have 3 atoms. Why do we see just one atom
per water molecule in this structure?**

-   we only see oxygen and not hydrogen due to the xray resolution is
    not depicted within the 3D digram for
    [1HSG](https://www.rcsb.org/structure/1hsg). When visualizing a
    molecule from a MOL file, the software used to display the structure
    may not render all of the atoms or bonds, depending on the settings
    or limitations of the software.

![](1HSG.png)

### pink sequences interact with the molecules usng the surroundings…

**Q5: There is a critical “conserved” water molecule in the binding
site. Can you identify this water molecule? What residue number does
this water molecule have**

-   yes this residue is idenfied HOH - 308 interacting with the two
    Asp’s at two amino acid interation at A-Asp 25, and B-Asp 25

**Q6: Generate and save a figure clearly showing the two distinct chains
of HIV-protease along with the ligand. You might also consider showing
the catalytic residues ASP 25 in each chain and the critical water (we
recommend *“Ball & Stick”* for these side-chains). Add this figure to
your Quarto document.**

**Discussion Topic:** Can you think of a way in which indinavir, or even
larger ligands and substrates, could enter the binding site?

-   Indinavir is a protease inhibitor used to treat HIV. It works by
    binding to the active site of HIV-1 protease, an enzyme required for
    the virus to replicate. The protease has a binding site that
    accommodates the substrate or ligand, such as indinavir. The process
    of ligand entry into the binding site is essential for the
    interaction between the enzyme and its inhibitor or substrate.

![](images/1HSG-7.png)

# 3. Introduction to Bio3D in R

<u>***install Bio3D***</u>:
[`Bio3D`](http://thegrantlab.org/bio3d/index.php) is an R package for
structural bioinformatics. Features include the ability to read, write
and analyze biomolecular structure, sequence and dynamic trajectory
data.

``` r
# Install Bio3D
# Load the Bio3D library
library(bio3d)
```

## Reading PDB file data into R

``` r
# Load a PDB file
pdb <- read.pdb("1hsg")
```

      Note: Accessing on-line PDB file

``` r
# Print information about the PDB file
pdb
```


     Call:  read.pdb(file = "1hsg")

       Total Models#: 1
         Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)

         Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
         Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)

         Non-protein/nucleic Atoms#: 172  (residues: 128)
         Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]

       Protein sequence:
          PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
          QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
          ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
          VNIIGRNLLTQIGCTLNF

    + attr: atom, xyz, seqres, helix, sheet,
            calpha, remark, call

**Q7: How many amino acid residues are there in this pdb object?**

-   There are 198 amino acid residues in this PDB object, as indicated
    by the line “Protein Atoms#: 1514 (residues/Calpha atoms#: 198)”.

**Q8: Name one of the two non-protein residues?**

-   One of the non-protein residues is HOH (water). This information can
    be found in the line “Non-protein/nucleic resid values: \[ HOH
    (127), MK1 (1) \]”.

**Q9: How many protein chains are in this structure?**

-   There are 2 protein chains in this structure, as shown by the line
    “Chains#: 2 (values: A B)”.

Note that the attributes `(+ attr:)` of this object are listed on the
last couple of lines. To find the attributes of any such object you can
use:

``` r
# Use attributes command to pull the last lines 
attributes(pdb)
```

    $names
    [1] "atom"   "xyz"    "seqres" "helix"  "sheet"  "calpha" "remark" "call"  

    $class
    [1] "pdb" "sse"

To access these individual attributes we use the `dollar-attribute`name
convention that is common with R list objects. For example, to access
the `atom` attribute or component use `pdb$atom`:

``` r
#use head command to pull frist 6 rows accessing the atom
head(pdb$atom)
```

      type eleno elety  alt resid chain resno insert      x      y     z o     b
    1 ATOM     1     N <NA>   PRO     A     1   <NA> 29.361 39.686 5.862 1 38.10
    2 ATOM     2    CA <NA>   PRO     A     1   <NA> 30.307 38.663 5.319 1 40.62
    3 ATOM     3     C <NA>   PRO     A     1   <NA> 29.760 38.071 4.022 1 42.64
    4 ATOM     4     O <NA>   PRO     A     1   <NA> 28.600 38.302 3.676 1 43.40
    5 ATOM     5    CB <NA>   PRO     A     1   <NA> 30.508 37.541 6.342 1 37.87
    6 ATOM     6    CG <NA>   PRO     A     1   <NA> 29.296 37.591 7.162 1 38.40
      segid elesy charge
    1  <NA>     N   <NA>
    2  <NA>     C   <NA>
    3  <NA>     C   <NA>
    4  <NA>     O   <NA>
    5  <NA>     C   <NA>
    6  <NA>     C   <NA>

## Predicting functional motions of a single structure

Let’s read a new PDB structure of Adenylate Kinase and perform Normal
mode analysis.

``` r
# Read the PDB file and store the structure in 'adk'
adk <- read.pdb("6s36")
```

      Note: Accessing on-line PDB file
       PDB has ALT records, taking A only, rm.alt=TRUE

``` r
# Print information about the loaded structure
adk
```


     Call:  read.pdb(file = "6s36")

       Total Models#: 1
         Total Atoms#: 1898,  XYZs#: 5694  Chains#: 1  (values: A)

         Protein Atoms#: 1654  (residues/Calpha atoms#: 214)
         Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)

         Non-protein/nucleic Atoms#: 244  (residues: 244)
         Non-protein/nucleic resid values: [ CL (3), HOH (238), MG (2), NA (1) ]

       Protein sequence:
          MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVT
          DELVIALVKERIAQEDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFDVPDELIVDKI
          VGRRVHAPSGRVYHVKFNPPKVEGKDDVTGEELTTRKDDQEETVRKRLVEYHQMTAPLIG
          YYSKEAEAGNTKYAKVDGTKPVAEVRADLEKILG

    + attr: atom, xyz, seqres, helix, sheet,
            calpha, remark, call

``` r
# Perform flexiblity prediction
m <- nma(adk)
```

     Building Hessian...        Done in 0.012 seconds.
     Diagonalizing Hessian...   Done in 0.261 seconds.

``` r
# Determine the class of 'm'
class(m)
```

    [1] "VibrationalModes" "nma"             

``` r
# Plot 'm'
plot(m)
```

![](gg_files/figure-gfm/unnamed-chunk-20-1.png)

``` r
# Create a trajectory file called 'adk_m7.pdb'
mktrj(m, file="adk_m7.pdb")
```

![](images/Screenshot%202023-05-09%20at%203.42.08%20PM.png)
