<!-- This file is used to create README.md
Note that the README.md document may need updating to change
'\<0.001' to '<0.001'. 
-->

 
# BayesAT Package


Bayesian adaptive algorithm implements multiple-stage interim analysis for the single-arm trial. The package includes a data-generating function, a Bayesian hypothesis testing function, and a Bayesian adaptive algorithm. 


## Description
This project is designed to...

## Table of Contents
- [Installation](#installation)
- [Contact](#contact)
- [Usage](#usage)
  
## Installation

1. Installing from CRAN:
   ```
   install.packages('BayesAT')
   ```
2. Installing the updated version from [GitHub](https://github.com/) with:
   ```{r}
   devtools::install_github("AdamYZhong/BayesAT")
   ```

## Contact
Maintainer email: [yuan.zhong@uhn.ca](mailto:yuan.zhong@uhn.ca).


## Usage

There are three main functions in BayesAT package:

1. `Simulate_Enroll()` generates multiple streams of data sets with survival time, censoring status, and enrollment time.

2. `Bayes_test()` conducts a hypothesis test through the Bayesian survival model. 

3. `BayesAT()` conducts or simulates Bayesian adaptive trials through multiple-stage interim analysis.

Two online documents provide detailed examples [Data Simulation](https://github.com/AdamYZhong/BayesAT/blob/main/vignettes/SimulationEnroll.pdf) and [Bayesian Adaptive Algorithm](https://github.com/AdamYZhong/BayesAT/blob/main/vignettes/BayesianInference.pdf).
