# TSBB-vs-JMBB

Code and materials for the work:  
**Bayesian joint modeling for Beta-Binomial longitudinal and survival data: a comparison between simultaneous and two-stage estimation** Galán-Arcicollar, C., Alvares,D., Najera-Zuloaga, J., & Lee, D.-J. (2025). 
## Overview

This repository contains the code used to perform the statistical analysis described in our work. We propose a Bayesian joint specification approach to joint modeling including beta-binomial distribution to analyze the relationship between patient-reported outcomes (PROs) and survival outcomes.
It contains the R and Stan code for all relevant analyses in the work. Real data and settings from the first scenario simulation (based on real data), are excluded because the data cannot be made publicly available.

##  Main Features

- **Beta-Binomial longitudinal modeling**: Incorporates a flexible and realistic modeling approach for patient-reported outcomes (PROs), which are often bounded and overdispersed.
  
- **Bayesian joint speficication to joint modeling framework**: Uses a Bayesian approach to joint modeling. It overcomes two-stage approaches proposed in the literature
  
- **Comparison of modeling approaches**:
  - The proposed **Beta-Binomial + Bayesian joint specification joint model**.
  - A **two-stage approach for joint modeling** assuming a Beta-Binomial longitudinal component.
  
- **Simulation study**: Assesses and compares the performance and bias of the different approaches under various scenarios.
- **Dynamic predictions**. The use of Bayesian approach allows us to perform dynamic predictions for new subjects.
