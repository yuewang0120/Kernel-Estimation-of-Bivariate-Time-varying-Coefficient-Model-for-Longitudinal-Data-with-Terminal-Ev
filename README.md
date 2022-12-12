# Kernel Estimation of Bivariate Time-varying Coefficient Model for Longitudinal Data with Terminal Event

# Author Contributions Checklist Form 

## Data

### Abstract

The dataset contains the Medicare claims from End-Stage Renal Disease (ESRD) patients, organized by the United States Renal Data System (USRDS). Our analysis only concerns patients that had their first ESRD service in the year of 2007. According to the data use agreement, the data can not be made available to the public. Researchers who want to obtain access to the data should contact the USRDS for details on the procedure for requesting access to the data.

### Availability

According to the data use agreement, the data cannot be made public. Those who wish to get access to the dataset could visit [USRDS SAFs](https://www.usrds.org/for-researchers/standard-analysis-files/) and submit a request. For reproducibility, we have created a dataset `usrds.RData` that mimics the real USRDS data set used in the paper. The pseudo dataset is publicly available under `data` folder at the time of submission.

### Description

See `usrds_desc.txt` under `data` folder.

## Code

### Abstract

Researchers could use the R scripts to reproduce the results in the paper (including Figures 1-3, Tables 1-4 and Figures S1-S2 in the online supplementary material). However, since Figure 3, Table 3-4 and Figure S1 are generated with the original dataset, they are not completely reproducible. Readers can access a version of Figure 3, Table 3-4 and Figure S1 produced using the pseudo dataset under “pseudo_figures_and_tables” folder of the repo, which are similar-looking and generated using the exact same procedure (will elaborate in the next section) as the corresponding figures and tables in the paper.

## Instructions for Use

### 1. The Main Simulation with All Time-varying coefficients (Section 4)

A bandwidth of 0.66 is used for simulation, this can be validated by replicating the cross-validation procedure

```r
source('code/simu_cv.R')
```

Run simulation with 1000 replications with the selected bandwidth $h=0.66$

```r
source('code/simu.R')
```

This will generate intermediate results, stored in `simu_result.RData` under `code` folder. Once this is done, run

```r
source('code/figure1-2.R')
```

to generate Figure 1 and 2 under the same folder.

### 2. A Separate Simulation with $\beta_3=0.5$ (Section 4)

A bandwidth of 0.68 is used for simulation, this can be validated by replicating the cross-validation procedure

```r
source('code/simu_cv_beta3=0.5.R')
```

Run simulation with 1000 replications with the selected bandwidth $h=0.68$

```r
source('code/simu_beta3=0.5.R')
```

This will generate intermediate results, stored in `simu_result_beta3=0.5.RData` under `code` folder. Once this is done, run

```r
source('code/table1-2.R')
```

to generate Table 1 and 2 under the same folder.

### 3. ESRD Application with A Reduced Model (Section 5, with Pseudo Dataset)

Replicate cross-validation

```r
source('code/pseudo_cv.R')
```

Fit the model with the selected bandwidth

```r
source('code/pseudo_fit.R')
```

This will generate the intermediate result `pseudo_fit.RData` under `code` folder. Run

```r
source('code/figure3.R')
```

to generate Figure 3. Run

```r
source('code/table3-4.R')
```

to generate Table 3 and 4.

### 4. ESRD Application with A Full Model (Section C in the Online Supplementary Material, with Pseudo Dataset)

Replicate cross-validation

```r
source('code/pseudo_full_cv.R')
```

Run

```r
source('code/pseudo_full_fit.R')
```

to generate intermediate result `pseudo_full_fit.RData` under `code` folder. Run

```r
source('code/figureS1.R')
```

to generate Figure S1. 

### 5. Simulation for the Locally Weighted Pseudo Likelihood Approach (Section D in the Online Supplementary Material)

Run

```r
source('code/likelihood_simu.R')
```

to generate intermediate result `likelihood_simu_result.RData` under `code` folder. Run

```r
source('code/figureS2.R')
```

to generate Figure S2.
