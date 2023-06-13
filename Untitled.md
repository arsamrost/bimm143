Class 8: Mini Lab Cancer
================
Arsam Rostami

## Preparing the data

For this we can use the `read.csv()` function to read the CSV
(comma-separated values) file containing the data 

``` r
# Save your input data file into your Project directory using WisconsinCancer.csv

fna.data <- "WisconsinCancer.csv"

# Complete the following code to input the data and store as wisc.df and use row.names() command 

wisc.df <- read.csv(fna.data, row.names=1)
```

Examine your input data with the use of `id` and `diagnosis` column. we
can use the `View()` or `head ()` function in next columns

``` r
# using the wisc.df command pull the data up
head(wisc.df)
```

             diagnosis radius_mean texture_mean perimeter_mean area_mean
    842302           M       17.99        10.38         122.80    1001.0
    842517           M       20.57        17.77         132.90    1326.0
    84300903         M       19.69        21.25         130.00    1203.0
    84348301         M       11.42        20.38          77.58     386.1
    84358402         M       20.29        14.34         135.10    1297.0
    843786           M       12.45        15.70          82.57     477.1
             smoothness_mean compactness_mean concavity_mean concave.points_mean
    842302           0.11840          0.27760         0.3001             0.14710
    842517           0.08474          0.07864         0.0869             0.07017
    84300903         0.10960          0.15990         0.1974             0.12790
    84348301         0.14250          0.28390         0.2414             0.10520
    84358402         0.10030          0.13280         0.1980             0.10430
    843786           0.12780          0.17000         0.1578             0.08089
             symmetry_mean fractal_dimension_mean radius_se texture_se perimeter_se
    842302          0.2419                0.07871    1.0950     0.9053        8.589
    842517          0.1812                0.05667    0.5435     0.7339        3.398
    84300903        0.2069                0.05999    0.7456     0.7869        4.585
    84348301        0.2597                0.09744    0.4956     1.1560        3.445
    84358402        0.1809                0.05883    0.7572     0.7813        5.438
    843786          0.2087                0.07613    0.3345     0.8902        2.217
             area_se smoothness_se compactness_se concavity_se concave.points_se
    842302    153.40      0.006399        0.04904      0.05373           0.01587
    842517     74.08      0.005225        0.01308      0.01860           0.01340
    84300903   94.03      0.006150        0.04006      0.03832           0.02058
    84348301   27.23      0.009110        0.07458      0.05661           0.01867
    84358402   94.44      0.011490        0.02461      0.05688           0.01885
    843786     27.19      0.007510        0.03345      0.03672           0.01137
             symmetry_se fractal_dimension_se radius_worst texture_worst
    842302       0.03003             0.006193        25.38         17.33
    842517       0.01389             0.003532        24.99         23.41
    84300903     0.02250             0.004571        23.57         25.53
    84348301     0.05963             0.009208        14.91         26.50
    84358402     0.01756             0.005115        22.54         16.67
    843786       0.02165             0.005082        15.47         23.75
             perimeter_worst area_worst smoothness_worst compactness_worst
    842302            184.60     2019.0           0.1622            0.6656
    842517            158.80     1956.0           0.1238            0.1866
    84300903          152.50     1709.0           0.1444            0.4245
    84348301           98.87      567.7           0.2098            0.8663
    84358402          152.20     1575.0           0.1374            0.2050
    843786            103.40      741.6           0.1791            0.5249
             concavity_worst concave.points_worst symmetry_worst
    842302            0.7119               0.2654         0.4601
    842517            0.2416               0.1860         0.2750
    84300903          0.4504               0.2430         0.3613
    84348301          0.6869               0.2575         0.6638
    84358402          0.4000               0.1625         0.2364
    843786            0.5355               0.1741         0.3985
             fractal_dimension_worst
    842302                   0.11890
    842517                   0.08902
    84300903                 0.08758
    84348301                 0.17300
    84358402                 0.07678
    843786                   0.12440

The first column is the `wisc.df$diagnosis` to determine which cell
samples are malignant or benign

``` r
# We can use -1 here to remove the first column
wisc.data <- wisc.df[,-1]
```

setup a separate new vector called `diagnosis` from the original data
set

``` r
# Create diagnosis vector for later 
diagnosis <- wisc.df[,1]
```

## Exploratory data analysis

Explore the data you created before (`wisc.data` and `diagnosis`) to
answer the following questions:

**Q1**. How many observations are in this data set?

``` r
# Use nrow command to determine the number of Obs. from wisc.df

nrow(wisc.df)
```

    [1] 569

-   Using `nrow()` command we can determine that there is **569
    observations**

**Q2**. How many of the observations have a malignant diagnosis?

``` r
# Load the data (replace "wisc.df" with the name of your data frame)
# Use table command on wis.df$diagnosis to pull the malignant and benign 
table(wisc.df$diagnosis)
```


      B   M 
    357 212 

``` r
# Create a new variable "diagnosis" containing the diagnosis column from the data frame
diagnosis<- wisc.df$diagnosis

# Remove the diagnosis column from the original data frame
wisc.df<- wisc.df[,-1]
```

-   Based on the `table()` command using `wisc.df$diagnosis` the is
    **357 benign** and **212 malignant** cells

**Q3**. How many variables/features in the data are suffixed
with `_mean`?

``` r
# use the length command and grep to pull the "_mean" and identify the colnames 

length(grep("_mean", colnames(wisc.df)))
```

    [1] 10

-   There are **10 variables** in the data when using `_mean` as a set
    parameter

# 2. Principal Component Analysis

## Performing PCA

The next step in your analysis is to perform principal component
analysis (PCA) on `wisc.data`

Check for the mean and std for the `wisc.data`, and will be using
`colMeans()` and `apply()`

``` r
# Check column means and standard deviations
colMeans(wisc.data)
```

                radius_mean            texture_mean          perimeter_mean 
               1.412729e+01            1.928965e+01            9.196903e+01 
                  area_mean         smoothness_mean        compactness_mean 
               6.548891e+02            9.636028e-02            1.043410e-01 
             concavity_mean     concave.points_mean           symmetry_mean 
               8.879932e-02            4.891915e-02            1.811619e-01 
     fractal_dimension_mean               radius_se              texture_se 
               6.279761e-02            4.051721e-01            1.216853e+00 
               perimeter_se                 area_se           smoothness_se 
               2.866059e+00            4.033708e+01            7.040979e-03 
             compactness_se            concavity_se       concave.points_se 
               2.547814e-02            3.189372e-02            1.179614e-02 
                symmetry_se    fractal_dimension_se            radius_worst 
               2.054230e-02            3.794904e-03            1.626919e+01 
              texture_worst         perimeter_worst              area_worst 
               2.567722e+01            1.072612e+02            8.805831e+02 
           smoothness_worst       compactness_worst         concavity_worst 
               1.323686e-01            2.542650e-01            2.721885e-01 
       concave.points_worst          symmetry_worst fractal_dimension_worst 
               1.146062e-01            2.900756e-01            8.394582e-02 

``` r
# Use the apply function to calculate the standard deviation for each column in the data frame
apply(wisc.data,2,sd)
```

                radius_mean            texture_mean          perimeter_mean 
               3.524049e+00            4.301036e+00            2.429898e+01 
                  area_mean         smoothness_mean        compactness_mean 
               3.519141e+02            1.406413e-02            5.281276e-02 
             concavity_mean     concave.points_mean           symmetry_mean 
               7.971981e-02            3.880284e-02            2.741428e-02 
     fractal_dimension_mean               radius_se              texture_se 
               7.060363e-03            2.773127e-01            5.516484e-01 
               perimeter_se                 area_se           smoothness_se 
               2.021855e+00            4.549101e+01            3.002518e-03 
             compactness_se            concavity_se       concave.points_se 
               1.790818e-02            3.018606e-02            6.170285e-03 
                symmetry_se    fractal_dimension_se            radius_worst 
               8.266372e-03            2.646071e-03            4.833242e+00 
              texture_worst         perimeter_worst              area_worst 
               6.146258e+00            3.360254e+01            5.693570e+02 
           smoothness_worst       compactness_worst         concavity_worst 
               2.283243e-02            1.573365e-01            2.086243e-01 
       concave.points_worst          symmetry_worst fractal_dimension_worst 
               6.573234e-02            6.186747e-02            1.806127e-02 

Use `prcomp()` function on the `wisc.data`, scaling if appropriate
output model to `wisc.pr`.

``` r
# Perform PCA on wisc.data by completing the following code
wisc.pr <- prcomp( wisc.data, scale=T)
```

Inspect a summary of the results with the `summary()` function.

``` r
# Look at summary of results
summary(wisc.pr)
```

    Importance of components:
                              PC1    PC2     PC3     PC4     PC5     PC6     PC7
    Standard deviation     3.6444 2.3857 1.67867 1.40735 1.28403 1.09880 0.82172
    Proportion of Variance 0.4427 0.1897 0.09393 0.06602 0.05496 0.04025 0.02251
    Cumulative Proportion  0.4427 0.6324 0.72636 0.79239 0.84734 0.88759 0.91010
                               PC8    PC9    PC10   PC11    PC12    PC13    PC14
    Standard deviation     0.69037 0.6457 0.59219 0.5421 0.51104 0.49128 0.39624
    Proportion of Variance 0.01589 0.0139 0.01169 0.0098 0.00871 0.00805 0.00523
    Cumulative Proportion  0.92598 0.9399 0.95157 0.9614 0.97007 0.97812 0.98335
                              PC15    PC16    PC17    PC18    PC19    PC20   PC21
    Standard deviation     0.30681 0.28260 0.24372 0.22939 0.22244 0.17652 0.1731
    Proportion of Variance 0.00314 0.00266 0.00198 0.00175 0.00165 0.00104 0.0010
    Cumulative Proportion  0.98649 0.98915 0.99113 0.99288 0.99453 0.99557 0.9966
                              PC22    PC23   PC24    PC25    PC26    PC27    PC28
    Standard deviation     0.16565 0.15602 0.1344 0.12442 0.09043 0.08307 0.03987
    Proportion of Variance 0.00091 0.00081 0.0006 0.00052 0.00027 0.00023 0.00005
    Cumulative Proportion  0.99749 0.99830 0.9989 0.99942 0.99969 0.99992 0.99997
                              PC29    PC30
    Standard deviation     0.02736 0.01153
    Proportion of Variance 0.00002 0.00000
    Cumulative Proportion  1.00000 1.00000

<u>**NOTE**</u>: This **wisc.df** will give us the data of **30
principle components** proposing std., prop, and cum. variance.

**Q4**. From your results, what proportion of the original variance is
captured by the first principal components (PC1)?

-   the proportion of the original variance for PC1 is **0.4427**

**Q5**. How many principal components (PCs) are required to describe at
least 70% of the original variance in the data?

-   To reach the **70%** of the original variance needs to reach from
    **PC1 to PC3**. so three principle components PC3

**Q6**. How many principal components (PCs) are required to describe at
least 90% of the original variance in the data?

-   To reach the **90%** of the original variance needs to reach from
    **PC1 to PC9**. So nine principle components PC9

## Interpreting PCA results

Create a biplot of the `wisc.pr` using the `biplot()` function.

``` r
# Using bioplot command on the wis.pr
biplot(wisc.pr)
```

![](Untitled_files/figure-gfm/unnamed-chunk-11-1.png)

**Q7.** What stands out to you about this plot? Is it easy or difficult
to understand? Why?

-   When using `biplot()` on `wisc.pr` , we are given a difficult plot
    with PC2 vs PC1 with two factors colored red and black which seems
    abstract trying to determine the benign and malignant cells
    initially observing

So lets generate a more standard scatter plot of each observation along
principal components 1 and 2 (i.e. a plot of PC1 vs PC2 available as the
first two columns of `wisc.pr$x`) and color the points by the diagnosis
(available in the `diagnosis` vector you created earlier).

``` r
# Scatter plot observations by components 1 and 2
plot( wisc.pr$x[,1] ,wisc.pr$x[,2], col= as.factor(diagnosis), xlab = "PC1", ylab = "PC2")
```

![](Untitled_files/figure-gfm/unnamed-chunk-12-1.png)

**Q8.** Generate a similar plot for principal components 1 and 3. What
do you notice about these plots?

``` r
# Scatter plot observations by components 1 and 3
plot( wisc.pr$x[,1] ,wisc.pr$x[,3], col= as.factor(diagnosis), xlab = "PC1", ylab = "PC3")
```

![](Untitled_files/figure-gfm/unnamed-chunk-13-1.png)

we can use the **ggplot2** package to make a more fancy figure we will
also need to add our `diagnosis` vector as a column if we want to use it
for mapping to the plot color aesthetic

``` r
# Create a data.frame for ggplot
df <- as.data.frame(wisc.pr$x)
df$diagnosis <- diagnosis

# Load the ggplot2 package
library(ggplot2)

# Make a scatter plot colored by diagnosis
ggplot(df) + 
  aes(PC1, PC2, col=as.factor(diagnosis)) + geom_point()
```

![](Untitled_files/figure-gfm/unnamed-chunk-14-1.png)

## Variance explained

Calculate the variance of each principal component by squaring the sdev
component of `wisc.pr` (i.e. `wisc.pr$sdev^2`). Save the result as an
object called `pr.var`.

``` r
# Calculate variance of each component
pr.var <- wisc.pr$sdev^2
head(pr.var)
```

    [1] 13.281608  5.691355  2.817949  1.980640  1.648731  1.207357

Assign this to a variable called `pve` and create a plot of variance
explained for each principal component.

``` r
# Variance explained by each principal component: pve
pve <- pr.var/sum(pr.var)

# Plot variance explained for each principal component
plot(pve, xlab = "Principal Component", 
     ylab = "Proportion of Variance Explained", 
     ylim = c(0, 1), type = "o")
```

![](Untitled_files/figure-gfm/unnamed-chunk-16-1.png)

### Alternative method…

``` r
# Alternative scree plot of the same data, note data driven y-axis
barplot(pve, ylab = "Precent of Variance Explained",
     names.arg=paste0("PC",1:length(pve)), las=2, axes = FALSE)
axis(2, at=pve, labels=round(pve,2)*100 )
```

![](Untitled_files/figure-gfm/unnamed-chunk-17-1.png)

## Few CRAN packages that are helpful for PCA

``` r
#install.packages("factoextra")
```

``` r
## ggplot based graph
#install.packages("factoextra")
library(factoextra)
```

    Welcome! Want to learn more? See two factoextra-related books at https://goo.gl/ve3WBa

``` r
fviz_eig(wisc.pr, addlabels = TRUE)
```

![](Untitled_files/figure-gfm/unnamed-chunk-19-1.png)

## Communicating PCA results

explain the mapping from the original features to the principal
components. The principal components are naturally ordered from the most
variance explained to the least variance explained.

**Q9.** For the first principal component, what is the component of the
loading vector (i.e. `wisc.pr$rotation[,1]`) for the
feature `concave.points_mean`? This tells us how much this original
feature contributes to the first PC.

``` r
# Extract the loading of the "concave.points_mean" variable for the first principal component (PC1)
wisc.pr$rotation["concave.points_mean",1]
```

    [1] -0.2608538

-   The <u>**-0.2608538**</u> indicates how much the
    ‘concave.points_mean’ feature contributes to the first principal
    component

# 3. Hierarchical clustering

As part of the preparation for hierarchical clustering, the distance
between all pairs of observations are computed. Furthermore, there are
different ways to link clusters together, with single, complete, and
average being the most common "linkage methods".

First scale the `wisc.data` data and assign the result to `data.scaled`.

``` r
# Scale the wisc.data data using the "scale()" function
data.scaled <- scale(wisc.data)
```

Calculate the (Euclidean) distances between all pairs of observations in
the new scaled data set and assign the result to`data.dist`.

``` r
# Calculate the Euclidean distance between rows in the scaled data
data.dist <- dist(data.scaled)
```

Create a hierarchical clustering model using complete linkage. Manually
specify the method argument to hclust() and assign the results
to `wisc.hclust`.

``` r
# Perform hierarchical clustering using complete linkage method
wisc.hclust <- hclust(data.dist, method = "complete")
```

## Results of hierarchical clustering

Let's use the hierarchical clustering model you just created to
determine a height

**Q10.** Using the `plot()` and `abline()` functions, what is the height
at which the clustering model has 4 clusters?

``` r
# Plot the dendrogram for the hierarchical clustering results
plot(wisc.hclust)

# Add a horizontal line at height 19 with red color and dashed line type
abline(h=19, col="red", lty=2)
```

![](Untitled_files/figure-gfm/unnamed-chunk-24-1.png)

-   For cluster, we are essentially ‘cutting’ the **Dendrogram** at that
    height. The number of vertical lines the horizontal line crosses
    corresponds to the number of clusters at that particular
    dissimilarity level. This can help you to determine the optimal
    number of clusters, based on your specific needs or the specific
    requirements of your analysis. The set **h value at 19** allows us
    to observe 4 clusters

## Selecting number of clusters

Use `cutree()` to cut the tree so that it has 4 clusters. Assign the
output to the variable `wisc.hclust.clusters`.

``` r
# stteing the parameters using the wis.clust set at h=19 with the cutree command 
wisc.hclust.clusters <- cutree(wisc.hclust, h=19)
```

We can use the `table()` function to compare the cluster membership to
the actual diagnoses.

``` r
# Use table to determine the four clusters of how many benign (B) and malignant(M) cells are present 
table(wisc.hclust.clusters, diagnosis)
```

                        diagnosis
    wisc.hclust.clusters   B   M
                       1  12 165
                       2   2   5
                       3 343  40
                       4   0   2

Here we picked four clusters and see that cluster 1 largely corresponds
to malignant cells (with `diagnosis` values of 1) whilst cluster 3
largely corresponds to benign cells (with `diagnosis` values of 0).

## Using different methods

number of different *"methods"* we can use to combine points during the
hierarchical clustering procedure. These
include `"single"`, `"complete"`, `"average"` and (my
favorite) `"ward.D2"`.

**Q12.** Which method gives your favorite results for the
same `data.dist` dataset? Explain your reasoning.

``` r
wisc.hclust <- hclust(data.dist, method = "complete")
plot(wisc.hclust)
```

![](Untitled_files/figure-gfm/unnamed-chunk-27-1.png)

``` r
wisc.hclust.avg <- hclust(data.dist, method = "average")
plot(wisc.hclust.avg)
```

![](Untitled_files/figure-gfm/unnamed-chunk-27-2.png)

``` r
wisc.hclust.sing <- hclust(data.dist, method = "single")
plot(wisc.hclust.sing)
```

![](Untitled_files/figure-gfm/unnamed-chunk-27-3.png)

``` r
wisc.hclust.D2 <- hclust(data.dist, method = "ward.D2")
plot(wisc.hclust.D2)
```

![](Untitled_files/figure-gfm/unnamed-chunk-27-4.png)

-   **wisc.hclust** with `method="ward.D2"` is the favorite due this
    method minimizes the total within-cluster variance. At each step,
    the pair of clusters with the smallest increase in total
    within-cluster variance after merging is chosen to merge. The method
    is suitable when you want to create compact and spherical clusters,
    and it is less prone to the chaining effect

# 4. Combining methods

## Clustering on PCA results

Note the two main branches of or dendrogram indicating two main
clusters - maybe these are malignant and benign. Let's find out!

``` r
# Perform hierarchical clustering on the first 7 principal components using the Ward's method
wisc.pr.hclust <- hclust(dist(wisc.pr$x[,1:7]), method="ward.D2")

# Cut the dendrogram into 2 groups and assign each sample to a group
grps <- cutree(wisc.pr.hclust, k=2)

# Create a table with the counts of samples in each group
table(grps)
```

    grps
      1   2 
    216 353 

``` r
# Create a contingency table comparing the groups obtained from clustering and the diagnosis
table(grps, diagnosis)
```

        diagnosis
    grps   B   M
       1  28 188
       2 329  24

``` r
# Create a scatter plot of the first two principal components colored by the cluster assignments
plot(wisc.pr$x[,1:2], col=grps)
```

![](Untitled_files/figure-gfm/unnamed-chunk-30-1.png)

``` r
# Create a scatter plot of the first two principal components colored by the cluster assignments with as.factor of diagnosis
plot(wisc.pr$x[,1:2], col=as.factor(diagnosis))
```

![](Untitled_files/figure-gfm/unnamed-chunk-31-1.png)

`hclust` cluster 1 is mostly "M" and cluster 2 is mostly "B" as we saw
from the results of calling `table(grps, diagnosis)`

``` r
# Convert the cluster assignments to a factor
g <- as.factor(grps)

# Get the levels of the factor
levels(g)
```

    [1] "1" "2"

Reversing the order…

``` r
# Now reverse the order using relevel function
g <- relevel(g,2)

# Get the levels of the factor
levels(g)
```

    [1] "2" "1"

``` r
# Plot using our re-ordered factor 
plot(wisc.pr$x[,1:2], col=g)
```

![](Untitled_files/figure-gfm/unnamed-chunk-34-1.png)

``` r
# Use the distance along the first 7 PCs for clustering i.e. wisc.pr$x[, 1:7]
wisc.pr.hclust <- hclust(dist(wisc.pr$x[,1:7]), 
                         method="ward.D2")
```

Cut this hierarchical clustering model into 2 clusters and assign the
results to `wisc.pr.hclust.clusters`.

``` r
# Cut the dendrogram into 2 clusters and assign each sample to a cluster
wisc.pr.hclust.clusters <- cutree(wisc.pr.hclust, k=2)
```

<br>

**Q13**. How well does the newly created model with four clusters
separate out the two diagnoses?

``` r
# Compare to actual diagnoses
table(wisc.pr.hclust.clusters,diagnosis)
```

                           diagnosis
    wisc.pr.hclust.clusters   B   M
                          1  28 188
                          2 329  24

**Q14**. How well do the hierarchical clustering models you created in
previous sections (i.e. before PCA) do in terms of separating the
diagnoses? Again, use the `table()` function to compare the output of
each model (`wisc.km$cluster` and `wisc.hclust.clusters`) with the
vector containing the actual diagnoses.

``` r
# Create a contingency table comparing the groups obtained from clustering (wisc.hclust.clusters) and the diagnosis
table(wisc.hclust.clusters, diagnosis)
```

                        diagnosis
    wisc.hclust.clusters   B   M
                       1  12 165
                       2   2   5
                       3 343  40
                       4   0   2

``` r
# Create a contingency table comparing the groups obtained from clustering (wis.pr.hclust.clusters) and the diagnosis
wis.pr.hclust.clusters <- cutree(wisc.pr.hclust, k=2)
table(wis.pr.hclust.clusters, diagnosis)
```

                          diagnosis
    wis.pr.hclust.clusters   B   M
                         1  28 188
                         2 329  24

``` r
# wis.pr.hclust.clusters vector analysis perscion Benign
(329/569)*100
```

    [1] 57.82074

``` r
# wis.pr.hclust.clusters vector analysis perscion Malignancy
(188/569)*100
```

    [1] 33.04042

-   57% accuracy for B and 33% accuracy for M collectively 90%

``` r
# wisc.hclust.clusters vector analysis perscion Benign
(343/569)*100
```

    [1] 60.2812

``` r
# wisc.hclust.clusters vector analysis perscion Malignancy
(165/569)*100
```

    [1] 28.99824

-   60% accuracy for B and 28% accuracy for M collectively 88%

-   Based on this we can conclude that wis.pr.hclust.clusters has the
    highest accuracy rate at 90% when both clusters are in summation but
    for treatment we want to lower the rate false-negative therefore
    having higher rate for Malignancy will help the treatment at 60% for
    wisc.hclust.clusters

# 6. Prediction

**Q16.** Which of these new patients should we prioritize for follow up
based on your results?

  
We will use the `predict()` function that will take our PCA model from
before and [new cancer cell
data](https://tinyurl.com/new-samples-CSV) and project that data onto
our PCA space.

``` r
#url <- "new_samples.csv"
# Read new samples data from the provided URL
url <- "https://tinyurl.com/new-samples-CSV"
new <- read.csv(url)

# Predict the principal components for the new samples using the PCA model
npc <- predict(wisc.pr, newdata=new)

# Print the predicted principal components
npc
```

               PC1       PC2        PC3        PC4       PC5        PC6        PC7
    [1,]  2.576616 -3.135913  1.3990492 -0.7631950  2.781648 -0.8150185 -0.3959098
    [2,] -4.754928 -3.009033 -0.1660946 -0.6052952 -1.140698 -1.2189945  0.8193031
                PC8       PC9       PC10      PC11      PC12      PC13     PC14
    [1,] -0.2307350 0.1029569 -0.9272861 0.3411457  0.375921 0.1610764 1.187882
    [2,] -0.3307423 0.5281896 -0.4855301 0.7173233 -1.185917 0.5893856 0.303029
              PC15       PC16        PC17        PC18        PC19       PC20
    [1,] 0.3216974 -0.1743616 -0.07875393 -0.11207028 -0.08802955 -0.2495216
    [2,] 0.1299153  0.1448061 -0.40509706  0.06565549  0.25591230 -0.4289500
               PC21       PC22       PC23       PC24        PC25         PC26
    [1,]  0.1228233 0.09358453 0.08347651  0.1223396  0.02124121  0.078884581
    [2,] -0.1224776 0.01732146 0.06316631 -0.2338618 -0.20755948 -0.009833238
                 PC27        PC28         PC29         PC30
    [1,]  0.220199544 -0.02946023 -0.015620933  0.005269029
    [2,] -0.001134152  0.09638361  0.002795349 -0.019015820

``` r
# Create a scatter plot of the first two principal components colored by the diagnosis
plot(wisc.pr$x[,1:2], col=as.factor(diagnosis))

# Add the new samples' predicted principal components to the plot with blue color, circle symbol, and larger size
points(npc[,1], npc[,2], col="blue", pch=16, cex=3)

# Add text labels to the new samples' points with white color
text(npc[,1], npc[,2], c(1,2), col="white")
```

![](Untitled_files/figure-gfm/unnamed-chunk-43-1.png)

patient 2 principle component acquires malignant cells along the PC1 and
PC2 due as indicated `as.factor(grps)` determining 1 as color black for
benign and 2 as color red for malignant cells. Hence patient 2 should be
more cautious due to encountering their cell within malignant cell as a
correlative analysis
