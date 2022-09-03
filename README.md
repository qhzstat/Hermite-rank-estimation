# Hermite-rank-estimation

This repo is for generating all the simulation outcomes and data examples in Section 5 in "A Method for Estimating the Hermite Rank of Long-memory Time Series."

## Data

1) Dow Jones Index (DJI) from $Yahoo!$ Finance and the Standard & Poor's 500 (S&P 500) from the Wall Street Journal, from 01 Jan 2010 to 01 Jan 2021 (cf. Stock_index_analysis.R).

3) Tree ring data from the International Tree-Ring Data Bank **[ITRDB](https://www.ncei.noaa.gov/products/paleoclimatology/tree-ring)** collected from Africa, Asia, Australia, Canada, Europe, Mexico, South America, and the USA, stored in the Standard Chronology File (*.crn) format (cf. Tree_ring_analysis.R).

## Code

- cache: this folder includes some simulation outcomes in *.rds format.
- data: this folder contains the S&P 500 data in SnP500_WSJ.csv and a crn_data folder for all the tree rings data.
- functions: Hermite_rank_estimator.R for the algo to estimate Hermite rank, LRD_Generator_new.R for generating Farima/Fractional Gaussin series, parallel.R for parallel computing.
- main: the folder contains all the code file to generate all the outcomes for simualtions and data examples seperately. 
  1) Tree_ring_analysis.R and Stock_index_analysis.R are for the two data examples, respectively. 
  2) Hermite_rank_simulations.R is for generating Table 1 & 2. 
  3) Hermite_rank_test_with_estimation.R is for generating Figure 6.
  4) Histogram_of_xbar.R is for generating Figure 1.
  5) JCGS_test.R is for generating Figure 3.
  6) phat_confidence_interval.R is for generating Figure 10 & 12.
- plots: this folder saves all the figures for the paper.
