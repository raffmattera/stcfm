# Replication package for "Forecasting house prices growth rates with factor models and spatio-temporal clustering"

This is a file explaining how to replicate the results of the paper "Forecasting house prices growth rates with factor models and spatio-temporal clustering".

Here you can find two main R codes:

- IJF_MF_simulations, needed for replicating the results (tables and figures) of the simulation study in Section 4;
- IJF_MF_empirical, needed for replicating the results of the forecasting experiment (Section 5).

The function implementing our methodology is called "STfacm" and requires the packages "TSclust" and "ClustGeo" to be used. The Ando and Bai (2017) approach is adapted to our setting from the "PDMIF" package.

# Simulations (Section 4):

To make the simulations work, shapefiles on US geography are needed. These are the files ."US_States_Boundaries" with different extensions. Ensure you have all the files with all the extensions in the folder where the R code is.

The figures in the R script will appear exactly as in the paper. About the tables in Section 4, the results will all be obtained from the R code, but we arranged the tables manually in the paper. Readers can verify that the numbers obtained with simulations (ARI, Forecast accuracy) correspond to those of the tables.

# Empirical (Section 5):

To replicate section 5, the following data files are needed:

- PriceUS.csv for the data on the house prices index;
- LatLongUS.xlsx for latitude and longitude (needed for clustering task in our algorithm);
- CPI_growth, LongTerm_ir and RealInc_growth are files including observable factors z.

  
