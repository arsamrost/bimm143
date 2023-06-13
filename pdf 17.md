class 17
================
Class 17 Mini-Project COVID-19 Vaccination Rates
05/31/2023

## Background

We will start by downloading the most recently dated *"Statewide
COVID-19 Vaccines Administered by ZIP Code"* CSV file from:

<https://data.ca.gov/dataset/covid-19-vaccine-progress-dashboard-data-by-zip-code>

# Getting Started

Be sure to move your downloaded CSV file to your project directory and
then read/import into an R object called `vax`. 

``` r
# Import vaccination data
vax <- read.csv("covid19vaccinesbyzipcode_test.csv")

# Show the first few rows
head(vax)
```

      as_of_date zip_code_tabulation_area local_health_jurisdiction      county
    1 2021-01-05                    94579                   Alameda     Alameda
    2 2021-01-05                    93726                    Fresno      Fresno
    3 2021-01-05                    94305               Santa Clara Santa Clara
    4 2021-01-05                    93704                    Fresno      Fresno
    5 2021-01-05                    94403                 San Mateo   San Mateo
    6 2021-01-05                    93668                    Fresno      Fresno
      vaccine_equity_metric_quartile                 vem_source
    1                              3 Healthy Places Index Score
    2                              1 Healthy Places Index Score
    3                              4 Healthy Places Index Score
    4                              1 Healthy Places Index Score
    5                              4 Healthy Places Index Score
    6                              1    CDPH-Derived ZCTA Score
      age12_plus_population age5_plus_population tot_population
    1               19192.7                20872          21883
    2               33707.7                39067          42824
    3               15716.9                16015          16397
    4               24803.5                27701          29740
    5               37967.5                41530          44408
    6                1013.4                 1199           1219
      persons_fully_vaccinated persons_partially_vaccinated
    1                       NA                           NA
    2                       NA                           NA
    3                       NA                           NA
    4                       NA                           NA
    5                       NA                           NA
    6                       NA                           NA
      percent_of_population_fully_vaccinated
    1                                     NA
    2                                     NA
    3                                     NA
    4                                     NA
    5                                     NA
    6                                     NA
      percent_of_population_partially_vaccinated
    1                                         NA
    2                                         NA
    3                                         NA
    4                                         NA
    5                                         NA
    6                                         NA
      percent_of_population_with_1_plus_dose booster_recip_count
    1                                     NA                  NA
    2                                     NA                  NA
    3                                     NA                  NA
    4                                     NA                  NA
    5                                     NA                  NA
    6                                     NA                  NA
      bivalent_dose_recip_count eligible_recipient_count
    1                        NA                        4
    2                        NA                        2
    3                        NA                        8
    4                        NA                        5
    5                        NA                        7
    6                        NA                        0
      eligible_bivalent_recipient_count
    1                                 4
    2                                 2
    3                                 8
    4                                 5
    5                                 7
    6                                 0
                                                                   redacted
    1 Information redacted in accordance with CA state privacy requirements
    2 Information redacted in accordance with CA state privacy requirements
    3 Information redacted in accordance with CA state privacy requirements
    4 Information redacted in accordance with CA state privacy requirements
    5 Information redacted in accordance with CA state privacy requirements
    6 Information redacted in accordance with CA state privacy requirements

**Q1.** What column details the total number of people fully vaccinated?

``` r
# persons_fully_vaccinated
colnames(vax[11])
```

    [1] "persons_partially_vaccinated"

-   **persons_fully_vaccinated** is the column details for people fully
    vaccinated

**Q2.** What column details the Zip code tabulation area?

``` r
# zip_code_tabulation_area
colnames(vax[2])
```

    [1] "zip_code_tabulation_area"

-   **zip_code_tabulation_area** is the column details for tabulation
    area

**Q3.** What is the earliest date in this dataset?

``` r
# convert the date column to Date class if it's not
vax$date <- as.Date(vax$as_of_date, format = "%Y-%m-%d")

# find the minimum (earliest) date
earliest_date <- min(vax$date, na.rm = TRUE)

# Print date
earliest_date
```

    [1] "2021-01-05"

-   The earliest date is **2021-01-05**

**Q4.** What is the latest date in this dataset?

``` r
# Find the maximum (latest) date
latest_date <- max(vax$date, na.rm = TRUE)

# Print date
latest_date
```

    [1] "2023-05-23"

-   The latest date is **2023-05-23**

As we have done previously, let's call the `skim()` function from
the **skimr** package to get a quick overview of this dataset:

``` r
# Apply the skim() function to your data
skimr::skim(vax)
```

|                                                  |        |
|:-------------------------------------------------|:-------|
| Name                                             | vax    |
| Number of rows                                   | 220500 |
| Number of columns                                | 20     |
| \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_   |        |
| Column type frequency:                           |        |
| character                                        | 5      |
| Date                                             | 1      |
| numeric                                          | 14     |
| \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_ |        |
| Group variables                                  | None   |

Data summary

**Variable type: character**

| skim_variable             | n_missing | complete_rate | min | max | empty | n_unique | whitespace |
|:--------------------------|----------:|--------------:|----:|----:|------:|---------:|-----------:|
| as_of_date                |         0 |             1 |  10 |  10 |     0 |      125 |          0 |
| local_health_jurisdiction |         0 |             1 |   0 |  15 |   625 |       62 |          0 |
| county                    |         0 |             1 |   0 |  15 |   625 |       59 |          0 |
| vem_source                |         0 |             1 |  15 |  26 |     0 |        3 |          0 |
| redacted                  |         0 |             1 |   2 |  69 |     0 |        2 |          0 |

**Variable type: Date**

| skim_variable | n_missing | complete_rate | min        | max        | median     | n_unique |
|:--------------|----------:|--------------:|:-----------|:-----------|:-----------|---------:|
| date          |         0 |             1 | 2021-01-05 | 2023-05-23 | 2022-03-15 |      125 |

**Variable type: numeric**

| skim_variable                              | n_missing | complete_rate |     mean |       sd |    p0 |      p25 |      p50 |      p75 |     p100 | hist  |
|:-------------------------------------------|----------:|--------------:|---------:|---------:|------:|---------:|---------:|---------:|---------:|:------|
| zip_code_tabulation_area                   |         0 |          1.00 | 93665.11 |  1817.38 | 90001 | 92257.75 | 93658.50 | 95380.50 |  97635.0 | ▃▅▅▇▁ |
| vaccine_equity_metric_quartile             |     10875 |          0.95 |     2.44 |     1.11 |     1 |     1.00 |     2.00 |     3.00 |      4.0 | ▇▇▁▇▇ |
| age12_plus_population                      |         0 |          1.00 | 18895.04 | 18993.87 |     0 |  1346.95 | 13685.10 | 31756.12 |  88556.7 | ▇▃▂▁▁ |
| age5_plus_population                       |         0 |          1.00 | 20875.24 | 21105.97 |     0 |  1460.50 | 15364.00 | 34877.00 | 101902.0 | ▇▃▂▁▁ |
| tot_population                             |     10750 |          0.95 | 23372.77 | 22628.50 |    12 |  2126.00 | 18714.00 | 38168.00 | 111165.0 | ▇▅▂▁▁ |
| persons_fully_vaccinated                   |     17711 |          0.92 | 14272.72 | 15264.17 |    11 |   954.00 |  8990.00 | 23782.00 |  87724.0 | ▇▃▁▁▁ |
| persons_partially_vaccinated               |     17711 |          0.92 |  1711.05 |  2071.56 |    11 |   164.00 |  1203.00 |  2550.00 |  42259.0 | ▇▁▁▁▁ |
| percent_of_population_fully_vaccinated     |     22579 |          0.90 |     0.58 |     0.25 |     0 |     0.44 |     0.62 |     0.75 |      1.0 | ▂▂▆▇▅ |
| percent_of_population_partially_vaccinated |     22579 |          0.90 |     0.08 |     0.09 |     0 |     0.05 |     0.06 |     0.08 |      1.0 | ▇▁▁▁▁ |
| percent_of_population_with_1\_plus_dose    |     23732 |          0.89 |     0.64 |     0.24 |     0 |     0.50 |     0.68 |     0.82 |      1.0 | ▂▂▅▇▆ |
| booster_recip_count                        |     74388 |          0.66 |  6373.43 |  7751.70 |    11 |   328.00 |  3097.00 | 10274.00 |  60022.0 | ▇▂▁▁▁ |
| bivalent_dose_recip_count                  |    159956 |          0.27 |  3407.91 |  4010.38 |    11 |   222.00 |  1832.00 |  5482.00 |  29484.0 | ▇▂▁▁▁ |
| eligible_recipient_count                   |         0 |          1.00 | 13120.40 | 15126.17 |     0 |   534.00 |  6663.00 | 22517.25 |  87437.0 | ▇▃▁▁▁ |
| eligible_bivalent_recipient_count          |         0 |          1.00 | 13016.51 | 15199.08 |     0 |   266.00 |  6562.00 | 22513.00 |  87437.0 | ▇▃▁▁▁ |

**Q5.** How many numeric columns are in this dataset?

``` r
# str(vax)
```

``` r
# Check which columns are numeric
numeric_columns <- sapply(vax, is.numeric)

# Count the number of numeric columns
num_numeric_columns <- sum(numeric_columns)

num_numeric_columns
```

    [1] 14

**Q6.** Note that there are "missing values" in the dataset. How
many `NA` values there in the`persons_fully_vaccinated` column?

``` r
# Count the number of NA values in the persons_fully_vaccinated column
num_na_values <- sum(is.na(vax$persons_fully_vaccinated))

# Print num_na_values
num_na_values
```

    [1] 17711

-   There is 17711 NA values within the
    `persons_fully_vaccinated` column

**Q7.** What percent of `persons_fully_vaccinated` values are missing
(to 2 significant figures)?

``` r
# Find the total number of values
total_values <- length(vax$persons_fully_vaccinated)

# Calculate the percentage of missing values
percentage_missing <- (num_na_values / total_values) * 100

# Print percentage_missing
percentage_missing
```

    [1] 8.0322

-   There is **8.03%** missing from the
    `persons_fully_vaccinated` values

**Q8.** \[Optional\]: Why might this data be missing?

-   many counties with their perspective zip codes have difficultly
    contacting and gathering data regarding COVID, therefore its creates
    difficulty to dteremine the accurate rate of vaccination in
    comparison to San Diego and LA

## Working with dates

One of the "character" columns of the data is `as_of_date`, which
contains dates in the Year-Month-Day format.

``` r
# Install lubridate if it's not already installed
# if (!require(lubridate)) {
# install.packages("lubridate")}

# Load the lubridate package
library(lubridate)
```


    Attaching package: 'lubridate'

    The following objects are masked from 'package:base':

        date, intersect, setdiff, union

What is today's date (at the time I am writing this obviously)

``` r
# Get today's date
today()
```

    [1] "2023-06-13"

-   today/s date **2023-05-31**

The `as_of_date` column of our data is currently not that usable. For
example we can't easily do math with it like answering the simple
question how many days have passed since data was first recorded:

``` r
# This will give an Error!
# today() - vax$as_of_date[1]
```

However if we convert our date data into a lubridate format things like
this will be much easier as well as plotting time series data later on.

``` r
# Specify that we are using the year-month-day format
vax$as_of_date <- ymd(vax$as_of_date)
```

Now we can do math with dates. For example: How many days have passed
since the first vaccination reported in this dataset?

``` r
# Calculate the difference
today() - vax$as_of_date[1]
```

    Time difference of 889 days

Using the last and the first date value we can now determine how many
days the dataset span?

``` r
# Calculate the difference
vax$as_of_date[nrow(vax)] - vax$as_of_date[1]
```

    Time difference of 868 days

**Q9.** How many days have passed since the last update of the dataset? 

``` r
# Load the lubridate package
library(lubridate)


# Ensure as_of_date is in Date format
vax$as_of_date <- ymd(vax$as_of_date)

# Calculate the number of days since the last update
days_since_last_update <- today() - max(vax$as_of_date, na.rm = TRUE)

# Print days_since_last_update
days_since_last_update
```

    Time difference of 21 days

-   There is **8 days** passed since the last update

**Q10.** How many unique dates are in the dataset (i.e. how many
different dates are detailed)?

``` r
# Calculate the number of unique dates
num_unique_dates <- length(unique(vax$as_of_date))

# Print num_unique_dates
num_unique_dates
```

    [1] 125

-   **125 unique** dates in the dataset

# Working with ZIP codes

One of the numeric columns in the dataset
(namely `vax$zip_code_tabulation_area`) are actually ZIP codes

``` r
# Install zipcodeR if it's not already installed
# if (!require(zipcodeR)) {
# install.packages("zipcodeR")}

# Load the zipcodeR package
library(zipcodeR)
```

``` r
# Geocode a zip code
geocode_zip('92037')
```

    # A tibble: 1 × 3
      zipcode   lat   lng
      <chr>   <dbl> <dbl>
    1 92037    32.8 -117.

Calculate the distance between the centroids of any two ZIP codes in
miles

``` r
zip_distance('92037','92109')
```

      zipcode_a zipcode_b distance
    1     92037     92109     2.33

More usefully, we can pull census data about ZIP code areas (including
median household income etc.). For example:

``` r
reverse_zipcode(c('92037', "92109") )
```

    # A tibble: 2 × 24
      zipcode zipcode_type major_city post_office_city common_city_list county state
      <chr>   <chr>        <chr>      <chr>                      <blob> <chr>  <chr>
    1 92037   Standard     La Jolla   La Jolla, CA           <raw 20 B> San D… CA   
    2 92109   Standard     San Diego  San Diego, CA          <raw 21 B> San D… CA   
    # ℹ 17 more variables: lat <dbl>, lng <dbl>, timezone <chr>,
    #   radius_in_miles <dbl>, area_code_list <blob>, population <int>,
    #   population_density <dbl>, land_area_in_sqmi <dbl>,
    #   water_area_in_sqmi <dbl>, housing_units <int>,
    #   occupied_housing_units <int>, median_home_value <int>,
    #   median_household_income <int>, bounds_west <dbl>, bounds_east <dbl>,
    #   bounds_north <dbl>, bounds_south <dbl>

# Focus on the San Diego area

We have two main choices on how to do this. The first using base R the
second using the **dplyr** package:

``` r
# Subset to San Diego county only areas
sd <- vax[vax$county == "San Diego", ]
```

Using **dplyr** the code would look like this:

``` r
# Load up dplyr library
library(dplyr)
```


    Attaching package: 'dplyr'

    The following objects are masked from 'package:stats':

        filter, lag

    The following objects are masked from 'package:base':

        intersect, setdiff, setequal, union

``` r
# Filter rows where county is "San Diego"
sd <- filter(vax, county == "San Diego")

# Count the number of rows in the resulting dataframe
nrow(sd)
```

    [1] 13375

Using **dplyr** is often more convenient when we are subsetting across
multiple criteria - for example all San Diego county areas with a
population of over 10,000.

``` r
# Filter rows where county is "San Diego" and age5_plus_population is greater than 10,000
sd.10 <- filter(vax, county == "San Diego" & age5_plus_population > 10000)
```

**Q11.** How many distinct zip codes are listed for San Diego County? 

``` r
# Number of distinct zip codes in San Diego County
num_distinct_zip_codes <- nrow(distinct(sd, zip_code_tabulation_area))

# print num_distinct_zip_codes
num_distinct_zip_codes
```

    [1] 107

**Q12.** What San Diego County Zip code area has the largest 12 +
Population in this dataset? 

``` r
#T idyverse
# Zip code area with the largest 12+ population
# largest_12plus_pop_zip = filter(sd, age5_plus_pop == max(sd$age5_plus_population))
# largest_12plus_pop_zip$

#Base R
largest_12plus_pop_zip = sd[sd$age5_plus_population == max (sd$age5_plus_population),]

# Use the unique command to determine zipcode
unique(largest_12plus_pop_zip$zip_code_tabulation)
```

    [1] 92154

-   The zip code of the county with the largest data is **92154**

``` r
# Zip code area with the largest 12+ population
largest_12plus_pop_zip <- sd %>%
  arrange(desc(age12_plus_population)) %>%
  head(n = 1) %>%
  pull(zip_code_tabulation_area)

# print largest_12plus_pop_zip
largest_12plus_pop_zip
```

    [1] 92154

-   Another way to approach and achieve the same result

**Q13.** What is the overall average "Percent of Population Fully
Vaccinated" value for all San Diego "County" as of "2023-05-31"? 

``` r
# Make sure the as_of_date column is in Date format
sd_may23 <- filter(sd, as_of_date == '2023-05-23')

# Calculate the average percent of population fully vaccinated
mean(sd_may23$percent_of_population_fully_vaccinated, na.rm = TRUE)
```

    [1] 0.7419654

-   There is **74.2%** of people are vaccinated

**Q14.** Using either ggplot or base R graphics make a summary figure
that shows the distribution of Percent of Population Fully Vaccinated
values as of "2023-02-28"?

``` r
# Load the ggplot2 package
library(ggplot2)

# Subset the data for 2023-05-31
sd_may23 <- filter(sd, county == 'San Diego' & 
                     as_of_date == '2023-05-23')
# Create a histogram
hist(sd_may23$percent_of_population_fully_vaccinated, 
     xlab = 'Percent of Population Fully Vaccinated', 
     main = 'Distribution of Percent of Population Fully Vaccinated
     San Diego, 2023-05-23',
     ylab =  "Count",)
```

![](pdf-17_files/figure-gfm/unnamed-chunk-31-1.png)

## Geom Option

``` r
# pull up the Library for the ggplot2
library(ggplot2)

# create a ggplot with geom_hoist in mind
ggplot(sd_may23) +
  aes(percent_of_population_fully_vaccinated) +
  geom_histogram(bins=12) +
  ggtitle('Distribution of Percent of Population Fully Vaccinated
     San Diego, 2023-05-23')
```

    Warning: Removed 8 rows containing non-finite values (`stat_bin()`).

![](pdf-17_files/figure-gfm/unnamed-chunk-32-1.png)

## Focus on UCSD/La Jolla

UC San Diego resides in the 92037 ZIP code area and is listed with an
age 5+ population size of 36,144.

``` r
# Filter rows for zip code 92037
ucsd <- filter(sd, zip_code_tabulation_area == "92037")

# Retrieve age5_plus_population value from the first row
age5_plus_pop <- ucsd[1, ]$age5_plus_population

age5_plus_pop
```

    [1] 36144

**Q15**. Using **ggplot** make a graph of the vaccination rate time
course for the 92037 ZIP code area:

### **Filtering UCSD By Zipcode**

``` r
# Create a ggplot using the function with using geom_line and _point as the parameters
ggplot(ucsd) +
  aes(as_of_date,
      percent_of_population_fully_vaccinated) +
  geom_point() +
  geom_line(group = 1) +
  ylim(c(0, 1)) +
  labs(title = "Vaccination Rate Time Course for ZIP Code 92037",
       x = "Date",
       y = "Percent Vaccinated")
```

![](pdf-17_files/figure-gfm/unnamed-chunk-34-1.png)

## Comparing to similar sized areas

Let's return to the full dataset and look across every zip code area
with a population at least as large as that of 92037
on*as_of_date* "2023-05-23".

``` r
# Subset to all CA areas with a population as large as 92037
vax.36 <- filter(vax, age5_plus_population > 36144 &
                as_of_date == "2023-05-23")

# To find the mean of 
mean(vax.36$percent_of_population_fully_vaccinated)
```

    [1] 0.7225892

**Q16**. Calculate the mean *"Percent of Population Fully
Vaccinated"* for ZIP code areas with a population as large as 92037 (La
Jolla) *as_of_date* "2023-02-28". Add this as a straight horizontal line
to your plot from above with the `geom_hline()` function?

``` r
# Calculate the mean percent vaccinated for ZIP codes with population as large as 92037 on 2023-05-23
mean_percent_vaccinated <- mean(vax.36$percent_of_population_fully_vaccinated, na.rm = TRUE)

# Plot the vaccination rate time course for ZIP code 92037
plot <- ggplot(ucsd) +
  aes(as_of_date, percent_of_population_fully_vaccinated) +
  geom_point() +
  geom_line(group = 1) +
  ylim(c(0, 1)) +
  labs(title = "Vaccination Rate Time Course for ZIP Code 92037",
       x = "Date",
       y = "Percent Vaccinated")

# Add a horizontal line for the mean percent vaccinated
plot + geom_hline(yintercept = mean_percent_vaccinated, linetype = "dashed", color = "green")
```

![](pdf-17_files/figure-gfm/unnamed-chunk-36-1.png)

**Q17.** What is the 6 number summary (Min, 1st Qu., Median, Mean, 3rd
Qu., and Max) of the *"Percent of Population Fully Vaccinated"* values
for ZIP code areas with a population as large as 92037 (La
Jolla) *as_of_date* "2023-02-28"?

``` r
# Calculate the 6-number summary
summary_stats <- summary(vax.36$percent_of_population_fully_vaccinated)

summary_stats
```

       Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
     0.3816  0.6469  0.7207  0.7226  0.7924  1.0000 

**Q18.** Using ggplot generate a histogram of this data.

``` r
# Load the ggplot2 package
library(ggplot2)

# Generate a histogram
ggplot(vax.36) +
  aes(percent_of_population_fully_vaccinated) +
  geom_histogram() +
  xlim(c(0,1))
```

    `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    Warning: Removed 2 rows containing missing values (`geom_bar()`).

![](pdf-17_files/figure-gfm/unnamed-chunk-38-1.png)

``` r
  labs(title = "Histogram of Percent of Population Fully Vaccinated in ZIP Code 92037",
       x = "Percent Vaccinated",
       y = "Count")
```

    $x
    [1] "Percent Vaccinated"

    $y
    [1] "Count"

    $title
    [1] "Histogram of Percent of Population Fully Vaccinated in ZIP Code 92037"

    attr(,"class")
    [1] "labels"

``` r
# genearte a histogram 
hist(vax.36$percent_of_population_fully_vaccinated, 
     xlab = 'Percent of Population Fully Vaccinated', 
     main = 'Distribution of Percent of Population Fully Vaccinated
     San Diego, 2023-05-23',
     ylab =  "Count",)
```

![](pdf-17_files/figure-gfm/unnamed-chunk-39-1.png)

**Q19**. Is the 92109 and 92040 ZIP code areas above or below the
average value you calculated for all these above?

``` r
vax %>% filter(as_of_date == "2023-02-28") %>%  
  filter(zip_code_tabulation_area=="92040") %>%
  select(percent_of_population_fully_vaccinated)
```

      percent_of_population_fully_vaccinated
    1                                0.55198

-   It is below the average and this code is using the tidyverse
    approach which is the pipes to do the filters. select is equivalent
    to selecting the variables. This is why Tidyverse is becoming very
    useful

**Q20.** Finally make a time course plot of vaccination progress for all
areas in the full dataset with a`age5_plus_population > 36144`.

``` r
# Load the ggplot2 package
library(ggplot2)

# Filter rows with age5_plus_population > 36144
vax.36.all <- filter(vax, age5_plus_population > 36144)

# Create the time course plot
ggplot(vax.36.all) +
  aes(as_of_date, percent_of_population_fully_vaccinated, group = zip_code_tabulation_area) +
  geom_line(alpha = 0.2, color = "purple") +
  ylim(0, 1) +
  labs(x = "Date", y = "Percent Vaccinated",
       title = "Vaccination Progress for Areas with Age 5+ Population > 36144",
       subtitle = "Time Course Plot") +
  geom_hline(yintercept = 0.722, linetype = "dashed", color = "red")
```

    Warning: Removed 185 rows containing missing values (`geom_line()`).

![](pdf-17_files/figure-gfm/unnamed-chunk-41-1.png)
