# grapsupOU-simulation-estimation-application
This repository contains R code to simulate and estimate graph supOU processes.

The code implements the new methodology introduced in the article "Statistical inference for Levy-driven graph supOU processes: From short- to long-memory in high-dimensional time series" by Shreya Mehta and Almut Veraart (Imperial College London).

In order to reproduce the results, please run the file GraphSupOU-Functions-Final.R first.
For the simulation study, the file Port-Sim.R needs to run first, followed by Port-Sim-Analyis.R.
The code for the empirical study is contained in EmpiricalStudy-Final.R and uses the RE Europe dataset introduced in Jensen, T., Pinson, P. RE-Europe, a large-scale dataset for modeling a highly renewable European electricity system. Sci Data 4, 170175 (2017). https://doi.org/10.1038/sdata.2017.175. Please download the file wind_signal_COSMO from that source first.
