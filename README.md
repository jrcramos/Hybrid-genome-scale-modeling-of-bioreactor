# Hybrid-genome-scale-modeling-of-bioreactor
Deep hybrid modeling of bioreactor cell culture data using neural networks embedded with a genome-scale model combined with first principles equations

**It will be a fork of https://github.com/jrcramos/Hybrid-modeling-of-bioreactor-with-LSTM with a genome scale model as spice!**

------------------------------------------------------------------------------------------------------------------------------------------------------

To run hybrid model training or simulation of two previously trained hybrid model
use the code "hybnet_train_main.m"
This codes shows how to define model structures and which parameter to use to run the
training process. Furthemore, several ways to plot and analyse the results

------------------------------------------------------------------------------------------------------------------------------------------------------
 This code "hybnet_train_main.m" comes with two predefined hybrid model structures one FFNN and
 one LSTM. Both structures were pre-trained and the data saved in  "hybrid_FFNN_1.mat"
 and "hybrid_LSTM_1.mat". It asks if the user wants to train these model structures again
 or simulate from saved files. Different plots are generated and a excel file
 "structures_fit_results.xlsx" with the overall results.

The folder ~/data contains data.xlsx with the feed DoE details and the simulations of 
concentrations over time generated using a dynamic model. This is a synthetic dataset,
the model was created based on the metabolic model proposed by Robitaille et al. (2015).
It also contains "data.mat" which is the import the concentrations, feed and other relevant 
informations. This file is generated using "/data/main_data_processing.m". When imported into 
matlab data(i).accum is the total amount of a metabolite that should be in the bioreactor over
time, it is the "sum of all added concentrations × volume added - sample volume × reactor concentration". 
The file also contains data(i).m_r which are the reacted amounts over time. 
The calculation of data(i).accum is made during the process of data(i).m_r calculation. 
The data(i).val are the estimated rates at different cell growth phases, it is estimated using "/data/RatesEstimation.m"
The latter is described in the supplementary material of this paper or the file "/data/read_me_data_processing.doc"
The folder ~/GEM contains the GEM of the synthetic ODE model used data generation and "main_least_square_regression.m" which
performes least square regression estimated rates vs predicted by the model) to show the model compatibility with the
estimated rates

------------------------------------------------------------------------------------------------------------------------------------------------------
