Class07: Machine Learning
================
Arsasm Rostami

## Example of K-means Clustering

First step, is to make up some data with the known structure, so we know
what the answer should be

``` r
# rnorm with distribution of zero and stnd dev. of 1
?rnorm
```

use the **`rnorm()`** function in R to achieve random numbers from a
normal distribution with mean zero and standard deviation of one

``` r
# Plot will give random values with std dev. of 1 with the mean of 0
plot(rnorm(10))
```

![](CLASS07_files/figure-gfm/unnamed-chunk-2-1.png)

``` r
# set the command to tmp to create a set of data with mean of -3,3
# Combine the random numbers and reverse the second set
tmp <- c(rnorm(30, mean = -3), rnorm(30, mean = 3))

# Plot the combined data
x <- cbind(x = tmp, y = rev(tmp))
plot(x)
```

![](CLASS07_files/figure-gfm/unnamed-chunk-3-1.png)

Now, we have some structured data in `x`. Let’s see if `k-means` is able
to identify the two groups

``` r
# per kmeans command to create an input for x, cneters, nstart
k <- kmeans(x, centers=2, nstart = 20)

# Print the k-means clustering result
print(k)
```

    K-means clustering with 2 clusters of sizes 30, 30

    Cluster means:
              x         y
    1  2.792546 -2.874796
    2 -2.874796  2.792546

    Clustering vector:
     [1] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1
    [39] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

    Within cluster sum of squares by cluster:
    [1] 61.15211 61.15211
     (between_SS / total_SS =  88.7 %)

    Available components:

    [1] "cluster"      "centers"      "totss"        "withinss"     "tot.withinss"
    [6] "betweenss"    "size"         "iter"         "ifault"      

Let’s explore `k`:

``` r
# use k with $ to inspect new elements for individual random data
# infinite random generation is required to get mean of 3
k$centers
```

              x         y
    1  2.792546 -2.874796
    2 -2.874796  2.792546

**`k$centers`** is a matrix where each row corresponds to a cluster
center, and the columns correspond to the variables (dimensions) of the
data. The number of rows is equal to the number of clusters specified in
the k-means algorithm

``` r
# use cluster command to color and Add cluster centroids to the plot
k$cluster
```

     [1] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1
    [39] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

``` r
# color the different groups with 2 different vectors, the commad plot will be able to seprate them
plot(x, col=k$cluster)
```

![](CLASS07_files/figure-gfm/unnamed-chunk-7-1.png)

Now, we can add the cluster centers using `col` to distinguish it from
other points and add `pch` package creates a box.

**NOTE**: Each `pch` package value <u>has different shapes</u>

``` r
# we need to first bring the cluster centers using points and setting colors
plot(x,col=k$cluster)

# Add cluster centroids to the plot
points(k$center, col = "blue", pch = 15)
```

![](CLASS07_files/figure-gfm/unnamed-chunk-8-1.png)

An Example when we select the wrong number of clusters

``` r
# Perform k-means clustering with 3 centers and 30 random starts
k_3 <-kmeans(x, centers = 3, nstart = 30)

# Plot the data with colors representing cluster membership
plot(x, col = k_3$cluster)
```

![](CLASS07_files/figure-gfm/unnamed-chunk-9-1.png)

## Example of Hierarchical Clustering

Let’s use the same data as before, which we stored in `x`. we will use
`hclust()` function

``` r
# this the input for the hclust using distance of x
clustering<-hclust(dist(x))

# print clustering to display the data showing the number of objects
clustering
```


    Call:
    hclust(d = dist(x))

    Cluster method   : complete 
    Distance         : euclidean 
    Number of objects: 60 

The algorithm has created two distinct groups based on the two elements
introduces

``` r
# advance command from the plot to create a hiearchical clustering
plot(clustering)
```

![](CLASS07_files/figure-gfm/unnamed-chunk-11-1.png)

Let’s add a horizontal line:

``` r
#over feeding the data for machine learning and not generalizing by over explaining
plot(clustering)

# Add horizontal line to demonstrate overfitting
abline(h=10,col= "pink")
```

![](CLASS07_files/figure-gfm/unnamed-chunk-12-1.png)

to get reults (i.e., membership vector) we need to cut the tree. The
function doing that is `cutree()`

``` r
# set subgroups eleement to cutree comand
subgroups <- cutree(clustering, h=10)

# print subgroups
subgroups
```

     [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2
    [39] 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2

Plotting this …

``` r
# Plot the data with colors representing subgroups
plot(x, col = subgroups)
```

![](CLASS07_files/figure-gfm/unnamed-chunk-14-1.png)

You can also “cut” tree with the number of clusters you want:

``` r
# define specificity using k groups
cutree(clustering, k=2)
```

     [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2
    [39] 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2

# Principle Component Analysis (PCA)

## \* PCA UK Food Data \*

Input the url file form `UK_foods.csv`

``` r
# paste the element url to download the data
url <- "https://tinyurl.com/UK-foods"
x <- read.csv(url)
head(x)
```

                   X England Wales Scotland N.Ireland
    1         Cheese     105   103      103        66
    2  Carcass_meat      245   227      242       267
    3    Other_meat      685   803      750       586
    4           Fish     147   160      122        93
    5 Fats_and_oils      193   235      184       209
    6         Sugars     156   175      147       139

**Q1**. How many rows and columns are in your new data frame named `x`?
What R functions could you use to answer this questions?

``` r
# Find the number of rows in the data frame x
num_rows <- nrow(x)
print(num_rows)
```

    [1] 17

``` r
# Find the number of columns in the data frame x
num_columns <- ncol(x)
print(num_columns)
```

    [1] 5

-   There is **17** number of rows and **5** number of columns

## Checking your data

For this task we can use the **View()** function to display all the data
(in a new tab in RStudio) or the**head()** and **tail()** functions to
print only a portion of the data

``` r
# Use View() fuction or Head() to display the data for the first 6 rows of x <chr>
head(x)
```

                   X England Wales Scotland N.Ireland
    1         Cheese     105   103      103        66
    2  Carcass_meat      245   227      242       267
    3    Other_meat      685   803      750       586
    4           Fish     147   160      122        93
    5 Fats_and_oils      193   235      184       209
    6         Sugars     156   175      147       139

Hmm, it looks like the row-names here were not set properly as we were
expecting 4 columns (one for each of the 4 countries of the UK - not 5
as reported from the **dim()** function).

``` r
# Set row names using the first column and remove the first column using minus indexing
# rownames(x) <- x[,1]
# x <- x[,-1]

# Print the modified data frame with row names and without the first column
# head(x)
```

This looks much better, now lets check the dimensions again:

``` r
# Define a function dim() that takes a matrix or data frame x as input
dim(x)
```

    [1] 17  5

**Much easier way to remove the unwanted the x \<chr\> from the
columns…**

``` r
#Lets change the parameters and use tow.names to perform same task 
x <- read.csv(url, row.names=1)
head(x)
```

                   England Wales Scotland N.Ireland
    Cheese             105   103      103        66
    Carcass_meat       245   227      242       267
    Other_meat         685   803      750       586
    Fish               147   160      122        93
    Fats_and_oils      193   235      184       209
    Sugars             156   175      147       139

> **Q2.** Which approach to solving the ‘row-names problem’ mentioned
> above do you prefer and why? Is one approach more robust than another
> under certain circumstances?
>
> -   Use a dedicated data structure for `row.names`: This approach
>     involves using a separate data structure to store the row names,
>     such as a named vector or a separate data frame. This allows the
>     row names to be accessed and manipulated independently of the data
>     itself. Using `row.names` is more robust comparing to other
>     methods, in addition the code is more flexible while comparing to
>     other method using a negative value such as `[ ,-1]` to remove a
>     vector. Also each time x \<- x\[,-1\] running this will remove a
>     vector which is not ideal. Specifically, the first run of the code
>     block would remove the first column of **`x`** and return a new
>     data frame with one less column. Subsequent runs of the code block
>     would remove the new first column of **`x`**, and so on, until all
>     columns of **`x`** have been removed.  

**Q3**: Changing what optional argument in the
above **barplot()** function results in the following plot?

-   when setting **`beside = FALSE`** produces the older behavior of
    stacking bars on top of each other. the **`beside`** argument
    controls whether bars should be drawn beside each other
    (**`beside = TRUE`**). To maintain backwards compatibility, the
    default behavior of **`barplot()`** was changed to draw bars beside
    each other, but setting **`beside = FALSE`**

Now we can generate some basic visualizations

``` r
# create a barplot for each country and use rainbow command to distinguish the 17 different catagories
barplot(as.matrix(x),col=rainbow(nrow(x)))
```

![](CLASS07_files/figure-gfm/unnamed-chunk-22-1.png)

Let’s refine the bar plot so it is more clear:

``` r
# Use beside command to distinguish the 17 different dimensions through the barplot command 
barplot(as.matrix(x),col=rainbow(nrow(x)), beside= TRUE)
```

![](CLASS07_files/figure-gfm/unnamed-chunk-23-1.png)

**Q5**: Generating all pairwise plots may help somewhat. Can you make
sense of the following code and resulting figure? What does it mean if a
given point lies on the diagonal for a given plot?

-   If a given point lies on the diagonal of a particular plot, it means
    that the x and y values of that point are equal. In the context of a
    scatterplot matrix, the diagonal plots would not provide any useful
    information, as they would simply be comparing a variable with
    itself. Instead, these diagonal plots are typically replaced with a
    histogram or density plot of the variable to display its
    distribution.

Other visualization that can be used…

``` r
# Use pairs() command which duplicates the plot allowing diagonal dimension
pairs(x, col= rainbow(nrow(x)), pch=16)
```

![](CLASS07_files/figure-gfm/unnamed-chunk-24-1.png)

**Q6**. What is the main differences between N. Ireland and the other
countries of the UK in terms of this data-set?

-   Using the diagonal dimension when comparing N. Ireland with other UK
    countries the clusters are not linear rather dispersed. Hence the
    `pair()` function is creating comparison and correlative analysis,
    the data suggests that the `pch` are not correlated with each other

## PCA to the rescue

Let’s apply PCA (principle command analysis). For that, we have to use
the command `prcomp()` . It is suspecting the Matrix transposed. This
function expects the transpose of our data (flipping the data).

``` r
# Transpose the data frame
t(x)
```

              Cheese Carcass_meat  Other_meat  Fish Fats_and_oils  Sugars
    England      105           245         685  147            193    156
    Wales        103           227         803  160            235    175
    Scotland     103           242         750  122            184    147
    N.Ireland     66           267         586   93            209    139
              Fresh_potatoes  Fresh_Veg  Other_Veg  Processed_potatoes 
    England               720        253        488                 198
    Wales                 874        265        570                 203
    Scotland              566        171        418                 220
    N.Ireland            1033        143        355                 187
              Processed_Veg  Fresh_fruit  Cereals  Beverages Soft_drinks 
    England              360         1102     1472        57         1374
    Wales                365         1137     1582        73         1256
    Scotland             337          957     1462        53         1572
    N.Ireland            334          674     1494        47         1506
              Alcoholic_drinks  Confectionery 
    England                 375             54
    Wales                   475             64
    Scotland                458             62
    N.Ireland               135             41

lets assigned this command to a proper name and use `prcomp()`

``` r
# transpose_matrix <- t(x)
# pca <- prcomp(transpose_matrix)
pca <- prcomp(t(x))
# depending on the object applied the summary command will give different principle components 
summary(pca)
```

    Importance of components:
                                PC1      PC2      PC3       PC4
    Standard deviation     324.1502 212.7478 73.87622 5.552e-14
    Proportion of Variance   0.6744   0.2905  0.03503 0.000e+00
    Cumulative Proportion    0.6744   0.9650  1.00000 1.000e+00

Let’s have a look by plotting the `pca` results

``` r
# create a plot for the pca command 
plot(pca)
```

![](CLASS07_files/figure-gfm/unnamed-chunk-27-1.png)

<u>***NOTE***</u>*: This graph is showing us something useful but the
components are not given*

We need to access the results of the PCA analysis

``` r
# Use atrributes command to define the characters using $names and $class
attributes (pca)
```

    $names
    [1] "sdev"     "rotation" "center"   "scale"    "x"       

    $class
    [1] "prcomp"

We can explore the data using pca\$x data frame:

``` r
#statistical technique used for data analysis to reduce the dimensionality of data while retaining as much information as possible.
pca$x
```

                     PC1         PC2         PC3           PC4
    England   -144.99315    2.532999 -105.768945  1.042460e-14
    Wales     -240.52915  224.646925   56.475555  9.556806e-13
    Scotland   -91.86934 -286.081786   44.415495 -1.257152e-12
    N.Ireland  477.39164   58.901862    4.877895  2.872787e-13

### Plotting:

``` r
# using [] and setting to ,1 and ,2 will set the parameters for the plot finding the components
plot(x=pca$x[,1],y=pca$x[,2] )
```

![](CLASS07_files/figure-gfm/unnamed-chunk-30-1.png)

**Q7**. Complete the code below to generate a plot of PC1 vs PC2. The
second line adds text labels over the data points.

``` r
# Add plot setting x to row 1 and y to row 2
plot(x=pca$x[,1],y=pca$x[,2], xlab="PC1", ylab="PC2", xlim=c(-270,500) )

# Add country names using colnames command and colors
text(x=pca$x[,1],y=pca$x[,2], colnames(x))
```

![](CLASS07_files/figure-gfm/unnamed-chunk-31-1.png)

**Q8.** Customize your plot so that the colors of the country names
match the colors in our UK and Ireland map and table at start of this
document.

``` r
# country names command using colname()
colnames(x)
```

    [1] "England"   "Wales"     "Scotland"  "N.Ireland"

``` r
# Add plot setting x to row 1 and y to row 2
plot(x=pca$x[,1],y=pca$x[,2], xlab="PC1", ylab="PC2", xlim=c(-270,500) )
# Setting the color parameters using arrow command to set them to different countries 
color_countries <- c("orange", "pink", "blue", "green")
# Add country names using colnames command and colors
text(x=pca$x[,1],y=pca$x[,2], colnames(x), col=color_countries)
```

![](CLASS07_files/figure-gfm/unnamed-chunk-33-1.png)

#### Below we can use the square of pca\$sdev , which stands for “standard deviation”, to calculate how much variation in the original data each PC accounts for

``` r
# Calculate the explained variance for each principal component
v <- round( pca$sdev^2/sum(pca$sdev^2) * 100 )

# Print the expalined variance 
v
```

    [1] 67 29  4  0

Or trying this method instead…

``` r
# Generate summary statistics for PCA
z <- summary(pca)

# Extract the importance of each principal component
z$importance
```

                                 PC1       PC2      PC3          PC4
    Standard deviation     324.15019 212.74780 73.87622 5.551558e-14
    Proportion of Variance   0.67444   0.29052  0.03503 0.000000e+00
    Cumulative Proportion    0.67444   0.96497  1.00000 1.000000e+00

The information can be summarized in a plot of the variances
(eigenvalues) with respect to the principal component number
(eigenvector number)

``` r
# Create a bar plot of the explained variance
barplot(v, xlab="Principal Component", ylab="Percent Variation")
```

![](CLASS07_files/figure-gfm/unnamed-chunk-36-1.png)

## Digging deeper (variable loadings)

``` r
# Lets focus on PC1 as it accounts for > 90% of variance 
par(mar=c(10, 3, 0.35, 0))
barplot( pca$rotation[,1], las=2 )
```

![](CLASS07_files/figure-gfm/unnamed-chunk-37-1.png)

**Q9**: Generate a similar ‘loadings plot’ for PC2. What two food groups
feature prominantely and what does PC2 maninly tell us about?

``` r
# Create a bar plot for loadings of PC2
par(mar = c(10, 3, 0.35, 0))
barplot(pca$rotation[, 2], las = 2, main = "Loadings for PC2", xlab = "Variables", ylab = "Loadings")
```

![](CLASS07_files/figure-gfm/unnamed-chunk-38-1.png)

-   Here we see observations (foods) with the largest positive loading
    scores that effectively “push” N. Ireland to right positive side of
    the plot `Fresh_potatoes`.

-   We can also see the observations/foods with high negative scores
    that push the other countries to the left side of the plot
    (including `Soft_drinks` ).

``` r
# The inbuilt biplot() can be useful for small datasets 
biplot(pca)
```

![](CLASS07_files/figure-gfm/unnamed-chunk-39-1.png)

## \*PCA of RNA-seq data\*

These data come from RNA-seq count data set to generate the data for the
expression of 2 conditions with WT in cell-lines and knockout. They are
defined as `wt` and `ko`

``` r
# download url for the dtataset to be loaded 
url2 <- "https://tinyurl.com/expression-CSV"
rna.data <- read.csv(url2, row.names=1)
head(rna.data)
```

           wt1 wt2  wt3  wt4 wt5 ko1 ko2 ko3 ko4 ko5
    gene1  439 458  408  429 420  90  88  86  90  93
    gene2  219 200  204  210 187 427 423 434 433 426
    gene3 1006 989 1030 1017 973 252 237 238 226 210
    gene4  783 792  829  856 760 849 856 835 885 894
    gene5  181 249  204  244 225 277 305 272 270 279
    gene6  460 502  491  491 493 612 594 577 618 638

**Q10**: How many genes and samples are in this data set?

``` r
# Define a function dim() that takes a matrix or data frame rna.data as input
dim(rna.data)
```

    [1] 100  10

-   I have 100 genes and 10 samples!

Let’s apply PCA:

``` r
# Perform PCA on RNA data
pca_rna= prcomp(t(rna.data), scale=TRUE)

# Display PCA summary
summary(pca_rna)
```

    Importance of components:
                              PC1    PC2     PC3     PC4     PC5     PC6     PC7
    Standard deviation     9.6237 1.5198 1.05787 1.05203 0.88062 0.82545 0.80111
    Proportion of Variance 0.9262 0.0231 0.01119 0.01107 0.00775 0.00681 0.00642
    Cumulative Proportion  0.9262 0.9493 0.96045 0.97152 0.97928 0.98609 0.99251
                               PC8     PC9      PC10
    Standard deviation     0.62065 0.60342 3.327e-15
    Proportion of Variance 0.00385 0.00364 0.000e+00
    Cumulative Proportion  0.99636 1.00000 1.000e+00

Let’s plot the principle command analysis 1 and 2

``` r
# Simple un polished plot of pc1 and pc2
plot(pca_rna$x[,1], pca_rna$x[,2], xlab = 'PC1', ylab = 'PC2')
```

![](CLASS07_files/figure-gfm/unnamed-chunk-43-1.png)

``` r
# Set up color samples
col_samples <- c(rep('blue',5), rep('red', 5))

#print the function
col_samples
```

     [1] "blue" "blue" "blue" "blue" "blue" "red"  "red"  "red"  "red"  "red" 

``` r
# Create a scatter plot of the first two principal components
plot(pca_rna$x[,1], pca_rna$x[,2], xlab = 'PC1', ylab = 'PC2', col= col_samples)
```

![](CLASS07_files/figure-gfm/unnamed-chunk-44-1.png)

A quick barplot summary of this Proportion of Variance for each PC can
be obtained by calling the `plot()` function directly on our prcomp
result object.

``` r
# Add a quick scree plot to demonstrate the variance witthin the rna
plot(pca_rna, main="Quick scree plot")
```

![](CLASS07_files/figure-gfm/unnamed-chunk-45-1.png)

``` r
## Variance captured per PC 
pca.var <- pca_rna$sdev^2

## Percent variance is often more informative to look at 
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca.var.per
```

     [1] 92.6  2.3  1.1  1.1  0.8  0.7  0.6  0.4  0.4  0.0

We can use this to generate our own scree-plot like this

``` r
# Assuming pca.var.per is an existing object containing the percent variation for each principal component
barplot(pca.var.per, main="Scree Plot", 
        names.arg = paste0("PC", 1:10),
        xlab="Principal Component", ylab="Percent Variation")
```

![](CLASS07_files/figure-gfm/unnamed-chunk-47-1.png)

### Using ggplot

We could use the **ggplot2** package here but we will first need a
data.frame as input for the main `ggplot()` function.

``` r
# Load required library
library(ggplot2)

# set the data as df for pca$x
df <- as.data.frame(pca_rna$x)

# Create a scatter plot of the first two principal components
ggplot(df) + 
  aes(PC1, PC2) + 
  geom_point()
```

![](CLASS07_files/figure-gfm/unnamed-chunk-48-1.png)

If we want to add a condition specific color and perhaps sample label
aesthetics for wild-type and knock-out samples we will need to have this
information added to our data.frame:

``` r
# Add a 'wt' and 'ko' "condition" column
df$samples <- colnames(rna.data) 
df$condition <- substr(colnames(rna.data),1,2)

# create another ggplot with legends as the geom label set to FALSE
p <- ggplot(df) + 
        aes(PC1, PC2, label=samples, col=condition) + 
        geom_label(show.legend = FALSE)

# print p
p
```

![](CLASS07_files/figure-gfm/unnamed-chunk-49-1.png)

And finally add some spit and polish

``` r
# Create a scatter plot of the first two principal components with specified title, subtitle, x, y labels, and caption
p + labs(title="PCA of RNASeq Data",
       subtitle = "PC1 clealy seperates wild-type from knock-out samples",
       x=paste0("PC1 (", pca.var.per[1], "%)"),
       y=paste0("PC2 (", pca.var.per[2], "%)"),
       caption="Class example data") +
     theme_bw()
```

![](CLASS07_files/figure-gfm/unnamed-chunk-50-1.png)

<u>**Note:**</u> By doing this analysis the change of expression can be
distinguished from each other. This method could be used for quality
control to detect the outliers using this technique

``` r
# Display the the barplot
barplot(pca_rna$rotation[,1])
```

![](CLASS07_files/figure-gfm/unnamed-chunk-51-1.png)

``` r
# use sort command for each gene
sort(pca_rna$rotation[,1])
```

          gene56       gene18        gene3       gene39       gene50       gene11 
    -0.103783479 -0.103774699 -0.103761385 -0.103744482 -0.103743341 -0.103719665 
          gene57       gene91        gene1       gene79       gene59       gene75 
    -0.103703675 -0.103698408 -0.103666005 -0.103639415 -0.103607438 -0.103592371 
          gene54       gene44       gene58       gene82       gene87       gene13 
    -0.103584153 -0.103504699 -0.103503980 -0.103481127 -0.103448562 -0.103399591 
          gene19       gene27       gene61       gene25       gene51       gene53 
    -0.103390599 -0.103374849 -0.103308945 -0.103302326 -0.103265591 -0.103245619 
          gene49       gene17       gene29       gene94       gene86       gene40 
    -0.103188532 -0.103013773 -0.102739689 -0.102692869 -0.102122719 -0.102003831 
          gene12       gene31       gene70       gene32       gene81       gene78 
    -0.102001924 -0.101768804 -0.101365212 -0.100677376 -0.100659777 -0.100499426 
          gene42       gene24       gene77       gene96       gene46       gene65 
    -0.098746675 -0.098284250 -0.097473626 -0.096658194 -0.096571619 -0.094219475 
          gene80       gene43       gene26        gene9       gene64       gene84 
    -0.093476477 -0.092001819 -0.085745836 -0.085460936 -0.078643996 -0.009263882 
          gene69       gene83        gene4        gene5       gene92       gene71 
     0.048197107  0.066065263  0.075320862  0.087428334  0.092534408  0.095664760 
          gene88        gene6       gene15       gene89       gene37        gene8 
     0.098226585  0.099670829  0.099993193  0.100038548  0.100467583  0.100759370 
          gene97       gene63       gene74       gene73       gene38       gene95 
     0.100787961  0.101468649  0.101747637  0.102001050  0.102080752  0.102142492 
          gene72       gene35       gene14       gene52       gene22       gene93 
     0.102347342  0.102382706  0.102478762  0.102519795  0.102725125  0.102950950 
          gene30       gene20       gene36       gene67       gene47       gene76 
     0.103044435  0.103121803  0.103412422  0.103453646  0.103502386  0.103514464 
           gene2       gene34       gene33       gene16        gene7       gene28 
     0.103514749  0.103525731  0.103592988  0.103598474  0.103609009  0.103638752 
          gene99       gene23       gene48       gene55       gene85       gene62 
     0.103649598  0.103681565  0.103682769  0.103695870  0.103698370  0.103713893 
          gene41       gene90       gene10       gene21       gene60       gene98 
     0.103716818  0.103777744  0.103783379  0.103787935  0.103805515  0.103837190 
          gene68       gene45       gene66      gene100 
     0.103839510  0.103840183  0.103845454  0.103870820 
