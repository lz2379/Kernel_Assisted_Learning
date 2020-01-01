This repository is for demostrating the kernel assisted learning (KAL) method proposed by Liangyu Zhu, Wenbin Lu, Michael R. Kosorko and Rui Song (2019).

#--------------------------------------------
File Explanations
#--------------------------------------------
'KAL-functions.R'  contains the functions for KAL method
'KAL-simulations.R' contains the code for applying KAL method to simulation settings 
'KAL-realdata.R' contains the code for applying KAL method to the IWPC dataset

The two methods below are for the purpose of comparison
'O-learning-functions.R'  contains the functions for O-learning method proposed by Guanhua Chen, Donglin Zeng, and Michael R. Kosorok (2016) (see paper for details)
'O-learning-simulations.R' contains the code for applying O-learning method to simulation settings 
'O-learning-realdata.R' contains the code for applying O-learning method to the IWPC dataset
'SVMW_0.2.tar.gz' is a package required for applying the O-learning method.

'discreteQ-functions.R'  contains the functions for discretized Q-learning method
'discreteQ-simulations.R' contains the code for applying discretized Q-learning method to simulation settings 
'discreteQ-realdata.R' contains the code for applying discretized Q-learning method to the IWPC dataset

'data-preprocess.R' contains the code for preprocessing the dataset
'Realdata-plots' contains the code for generating plots in the dataset

'/Figures/' contains the figures generated for the paper
'/results/' contains the results from simulation studies and real data application in 'csv' format.


#--------------------------------------------
IWPC Dataset
#--------------------------------------------
The original data: 'iwpc_data_7_3_09_revised.csv' is public available:  It can be downloaded from https://www.pharmgkb.org/downloads/ and it is the "IWPC Data"
The processed data: 'data-preprocessed.csv'

