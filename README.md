# Hybrid-modeling-of-bioreactor-with-LSTM
Deep hybrid modeling of a cell culture bioreactor data using Long Short-Term Memory (LSTM) networks combined with first principles equations

------------------------------------------------------------------------------------------------------------------------------------------------------

**To run hybrid model training or simulation of two previously trained hybrid model
use the code "hybnet_train_main.m"
This codes shows how to define model structures and which parameter to use to run the
training process. Furthemore, several ways to plot and analyse the results**

------------------------------------------------------------------------------------------------------------------------------------------------------
 This code "hybnet_train_main.m" comes with two predefined hybrid model structures one Feed Forward Neural Network (FFNN) and
 one Long Short-Term Memory (LSTM). Both structures were pre-trained and the data saved in  "hybrid_FFNN_1.mat"
 and "hybrid_LSTM_1.mat". It asks if the user wants to train these model structures again
 or simulate from saved files. Different plots are generated and a excel file
 "structures_fit_results.xlsx" with the overall results.

 The folder ~/data contains data.xlsx with the feed DoE details and the
 simulations of concentrations over time generated using a dynamic model. 
 This is a synthetic dataset, the model was created based on the metabolic model 
 proposed by Robitaille et al. (2015). It also contains "data.mat"
 which is the import the concentrations, feed and other relevant informations. 
 This file is generated using "/data/main_data_processing.m".
 When imported into matlab data(i).accum is the total
 amount of a metabolite that should be in the bioreactor over time, it is
 the "sum of all added concentrations × volume added - sample volume × reactor concentration". The
 file also contains data(i).m_r which are the reacted amounts over
 time. The calculation of data(i).accum is made during the process of
 data(i).m_r calculation. The latter is described in the supplementary
 material of this paper or the file "/data/read_me_data_processing.doc"


------------------------------------------------------------------------------------------------------------------------------------------------------
Data for the paper: 
Deep hybrid modeling of a HEK293 process: combining Long Short-Term Memory (LSTM) networks with first principles equations												
João R. C. Ramos, José Pinto, Gil Poiares-Oliveira, Ludovic Peeters, Patrick Dumas, Rui Oliveira

How to cite this article: 
Ramos, J. R. C., Pinto, J., Poiares‐Oliveira, G., Peeters, L., Dumas, P., & Oliveira, R. (2024). Deep hybrid modeling of 
a HEK293 process: Combining long short‐term memory networks with first principles equations.
Biotechnology and Bioengineering, 1–15.
https://doi.org/10.1002/bit.28668
-----------------------------------------------------------------------------------------------------------------------------------------------------
 
