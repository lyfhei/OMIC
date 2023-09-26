# Orthogonal Multimodality Integration and Clustering (OMIC)

## 1. Introduction 
OMIC is a method that apply projection method to integrate multiple sources of information while accounting for the dependence among them.
It excels at modeling the relationships among multiple variables, facilitating scalable computation,
and preserving accuracy in cell clustering compared to existing methods



## 2. Environment 
R version 4.0 or greater is required 
### Required packages: 
The detailed requirements can be found [here]([https://github.com/JiazhangCai/WEST/blob/main/requirements.txt](https://github.com/lyfhei/OMIC/blob/main/requirements)).

## 3.Tutorial 
The [tutorial](https://github.com/JiazhangCai/WEST/blob/main/tutorial.ipynb) provides a pipeline of how to implement
WEST on an actual sample, including the description for every parameter. 

## 4. Data used in the paper 
The DLPFC data used in the paper can be found [here](http://research.libd.org/spatialLIBD/).     
The mouse embryo data used in the paper can be found [here](https://crukci.shinyapps.io/SpatialMouseAtlas/).      
The simulated data used in the paper can be found in the folder [sim_data](https://github.com/JiazhangCai/WEST/tree/main/sim_data).       
The reference-based synthetic data is already in the [folder](https://github.com/JiazhangCai/WEST/tree/main/sim_data/ref_based).        
The reference-free synthetic data can be generated from [shiny.rdata](https://github.com/JiazhangCai/WEST/blob/main/sim_data/ref_free/shiny.rdata) using 
the code [simulation_generate.R](https://github.com/JiazhangCai/WEST/blob/main/sim_data/ref_free/shiny.rdata).
