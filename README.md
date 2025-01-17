# Replication files for "Forecasting house prices growth rates with factor models and spatio-temporal clustering"

Here you can find the files to replicate the results of the paper "Forecasting house prices growth rates with factor models and spatio-temporal clustering", written by Raffaele Mattera and Philip Hans Franses.

The results can be obtained from the following two main R codes (R version: 4.4.1):

- IJF_MF_simulations, for replicating the results (tables and figures) of the simulation study (Section 4);
- IJF_MF_empirical, to replicate the results of the forecasting experiment (Section 5).

The function implementing our methodology is called "STfacm" and requires the packages "TSclust" and "ClustGeo" to be used. The Ando and Bai (2017) approach is adapted from the "PDMIF" package. 

The code "IJF_MF_simulations.R" requires the following package versions:

ggplot2_3.5.1; maptools_0.9-5; sp_2.1-4; sf_1.0-16; ARDL_0.2.4; splm_1.6-5; simstudy_0.8.0;  TSclust_1.3.1;  cluster_2.1.6; pdc_1.0.3; ClustGeo_2.1;
forecast_8.23.0; sandwich_3.1-0;  lmtest_0.9-40; zoo_1.8-12; rio_1.1.1

!! We note that the package "maptools" has been archived from the CRAN on 2023-10-16. Users can either download and install the version "maptools_0.9-5" from the CRAN archive (here: https://cran.r-project.org/src/contrib/Archive/maptools/) or upload the "ushp.RDS" file included in this repository. The "ushp.RDS" object, which can be obtained via the maptools package as indicated in the code, includes information about the US geography needed for the simulation study. This RDS file can be imported using the "readRDS" function of the base package. In the "IJF_MF_simulations.R" code both solutions are available.

The code IJF_MF_empirical.R" requires the following package versions:

 ggplot2_3.5.1;  panelvar_0.5.5;  tidyr_1.3.1; dplyr_1.1.4;  ARDL_0.2.4; splm_1.6-5; dynlm_0.3-6; POET_2.0; TSclust_1.3.1; cluster_2.1.6; pdc_1.0.3; 
 ClustGeo_2.1; forecast_8.23.0; sandwich_3.1-0;  lmtest_0.9-40; zoo_1.8-12; rio_1.1.1

To reproduce the results correctly, we suggest downloading data and codes in the same folder (e.g. as a zip file) and running the R codes directly from this folder.

Please, feel free to contact the corresponding author (raffaele.mattera@uniroma1.it) for any questions on the paper and its code.

# Simulations (Section 4):

It is important to load the shapefiles for the US geography to obtain the simulation results. You can find the files called ."US_States_Boundaries" with different extensions. Please, ensure that all these files are in the same folder where the R code is executed.

G3_dgp1, G3_dgp2, G4_dgp1, G4_dgp2 are lists containing the results of all the 1000 simulations, considering 3 or 4 clusters under DGP1 or DGP2. Tables and figures in Section 4 are extracted from these lists. The results in Tables 1 and 2 will be obtained from the R code "IJF_MF_simulations", but we arranged the tables manually in the paper. Therefore, the readers can verify that the numbers obtained with simulations (ARI, Forecast accuracy) correspond to those of the paper, even if the presentation may be slightly different. About the Figures 2 and 3, we provided a code to obtain a subfigure with a generic time series length T. The readers are required to just change the length of the simulated time series T (called "tlength") in the code to obtain other subfigures.

# Empirical (Section 5):

To replicate section 5, the following data files are needed:

- PriceUS.csv for the data on the house prices index;
- LatLongUS.xlsx for latitude and longitude (needed for clustering task in our algorithm);
- CPI_growth, LongTerm_ir and RealInc_growth are files including observable factors z.

Figure 5 of the paper (showing the clustering results using full-sample data) is obtained outside the R environment, specifically from the website https://www.mapchart.net/usa.html and using the clustering results. The other figures and tables will be generated directly from the R code "IJF_MF_empirical". We indicate with comments when a given part of the code is used to generate a specific table or figure. We use the same numbering of tables and figures as appears in the manuscript.
