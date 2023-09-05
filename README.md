# Fuzzy Classification based on Depth Function (FC-DF)

 This repository contains the R code to apply the Fuzzy Classification based on Depth Function (FC-DF) algorithm. Besides, examples to ease the application of FC-DF are included:
 - example1.R contains a simple application of the FC-DF with given hyper-parameters values.
 - example2.R extends the first example including the search for suitable hyper-parameter values.
 
  In addition, public data sets used in the work  *Fuzzy classification with distance-based depth prototypes: High-dimensional Unsupervised and/or Supervised problems*, as well as the R files (AlizadehReproducibility.R, ClevelandReproducibility.R) including the performed analysis are included. A docker container is provided, and isntructions for its downloading and execution are given in a subsequent section.

The source of the public data set are:
 - A. A. Alizadeh, M. B. Eisen, R. E. Davis, et al., Distinct types of diffuse large b-cell lymphoma identified by gene expression profiling, *Nature*, (2000), 403, 503–511
 - D. Dua, C. Graff, UCI machine learning repository (2017). URL http://archive.ics.uci.edu/ml
 
## Running the examples
 
 You might need to install the libraries *dplyr*, *fclust*, *mnormt*, *ICGE*, *splitTools* if not already present in your R installation.

```bash
root@fb8b4292bc50:/root# R

R version 4.2.2 (2022-10-31) -- "Innocent and Trusting"
Copyright (C) 2022 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

  Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> install.packages(c("dplyr", "fclust", "mnormt", "ICGE", "splitTools"))
```

 
## Docker container
 
 
