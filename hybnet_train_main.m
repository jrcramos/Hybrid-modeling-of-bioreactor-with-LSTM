%---------------------------------------------------------------------------
% Deep hybrid modeling of a HEK293 process: combining Long Short-Term Memory (LSTM) networks with first principles equations												
% João R. C. Ramos1, #, José Pinto1, #, Gil Poiares-Oliveira1, Ludovic Peeters2, Patrick Dumas2, Rui Oliveira1																							
% 1 LAQV-REQUIMTE, Department of Chemistry, NOVA School of Science and Technology, NOVA University Lisbon, 2829-516 Caparica, Portugal 												
% 2 GSK, 89 rue de l'Institut, 1330 Rixensart, Belgium																							
% *Correspondence:												
% Rui Oliveira												
% rmo@fct.unl.pt
%---------------------------------------------------------------------------


%---------------------------------------------------------------------------
% This codes comes with two predefined hybrid model structures one FFNN and
% one LSTM. Both structures were pre-trained and the data saved in  "hybrid_FFNN_1.mat"
% and "hybrid_LSTM_1.mat". It asks if the user wants to train these model structures again
% or simulate from saved files. Different plots are generated and a excel file
% "structures_fit_results.xlsx" with the overall results.
%
% The folder ~/data contains data.xlsx with the feed DoE details and the
% simulations of concentrations over time generated using a dynamic model. 
% This is a synthetic dataset, the model was created based on the metabolic model 
% proposed by Robitaille et al. (2015). It also contains "data.mat"
% which import the concentrations, feed and other relevant informations.
% Futhermore when imported into matlab data(i).accum is the total
% amount of a metabolite that should be in the bioreactor over time, it is
% the "sum of all added concentrations × volume added - sample volume × reactor concentration". 
% The file also contains data(i).m_r which are the reacted amounts over
% time. The calculation of data(i).accum is made during the process of
% data(i).m_r calculation. The latter is described in the supplementary
% material of this paper
%---------------------------------------------------------------------------

%---------------------------------------------------------------------------
% DoE for experiment Br1 - Br9 (see data.xlsx):
%---------------------------------------------------------------------------
%
% 			       Br5		
%                   |
% 					|
% 		 Br1 -------------- Br2		
% 		  |					 |
% 		  |                  |
% Br7 ----|	       Br9		 |---- Br8
% 		  |					 |
% 		  |					 |
% 		 Br3 -------------- Br4		
% 					|			
% 					|			
% 			       Br6					
%


clc
close all
w = warning ('off','all');

%---------------------------------------------------------------------------
% define the layers of each model structures, names to save their results to, number of
% iterations to use (niter), number of repetitions (nruns), number of
% principal compoenents for PCA of reacted amounts (npcs), which
% experiments to use for training(Indtr)/validation(Indcr)/test(Indtt,optional)
%---------------------------------------------------------------------------
structures={};
niter=40000; % fixed for comparability
nruns=1; % DoE data so only one training repetition performed, in this isntance nruns may be changed,
         % but the same permutation of experiments [1 2 3 4 9] will be used
         % (repmat([1 2 3 4 9],nruns,1))for training, because of known DoE conditions
npcs= 4; % number of pricipal compenent for reaction correlation matrix
Indtr=repmat([1 2 3 4 9],nruns,1); % use ccdesign DoE square for hybnet training
Indcr=repmat([7],nruns,1); % use of diamond edge of ccdesign DoE for validation
Indtt=repmat([5 6 8],nruns,1); % use other 3 diamond edge for testing


%---------------------------------------------------------------------------
% simulate or run training? training is performed with the hybrid
% structures defined bellow
%---------------------------------------------------------------------------
list = {'Simulate','New training'};
[indx,tf] = listdlg('PromptString','Train hybrid NN or simulate:',...
                     'SelectionMode','single',...
                      'ListString',list);
if indx==2 && tf==1
    
    % hybrid FFNN
    structures{1,1}     = {hnetfflayer(27,10,'relu')  hnetfflayer(10,10,'relu') hnetfflayer(10,npcs,'relu') }; % FFNN hybrid structure
    structures{1,2}='hybrid_FFNN_1';% name of save file
    structures{1,3}=niter;structures{1,4}=nruns;structures{1,5}=npcs;
    
    % Hybrid LSTM
    structures{end+1,1}     = {hnetfflayer(27,10,'relu')  hnetfflayer(10,10,'relu') hnetLSTMlayer(10,npcs) }; % LSTM hybrid structure
    structures{end,2}='hybrid_LSTM_1';% name of save file
    structures{end,3}=niter;structures{end,4}=nruns;structures{end,5}=npcs;
    
    
    % starting training process, use "parfor" instead of "for" to run several trainings parallel 
    for i=1:size(structures,1)
%         close all
        hybnet_train_fun(structures{i,2},structures{i,1},Indtr,Indcr,structures(i,3:end))
        
    end
else
    hybnet_train_fun('hybrid_LSTM_1'); % simulate saved LSTM model and do plots
    hybnet_train_fun('hybrid_FFNN_1'); % simulate saved FFNN model and do plots
    
    %%%%%%%%%%%%%%%%%% other plots %%%%%%%%%%%%%%%%%%
    
    %%% plots BR1 and Br6 concentrations data vs model simulations
    plot_concentrations_paper('hybrid_LSTM_1',4,[1 6]); 
    
    %%% plots normalized predictions vs experimental concentrations for each 
    %%% process variable indivdually of Br1 to Br9 with R^2
    plot_concentrations_indiv_paper('hybrid_LSTM_1',4,1:9); 
    
    %%% plots the rates of Br1 and Br6
    plot_predicted_rates_paper('hybrid_LSTM_1',4,[1 6])
end


%---------------------------------------------------------------------------
% generate a table with the overall results and save it to an excel file
%---------------------------------------------------------------------------
if indx==2 && tf==1
    results=[];
    for i=1:size(structures,1)
        count=2;
        load([structures{i,2} '.mat']) % load saved results
        table_dnn={'hybrid structure' 'dropout rate' 'minibatch size' 'alpha'	'alphaf' 'Niter multiplier' 'RMSE_tr' 'RMSE_tt' ...
                   'RMSE_all' 'MSE_tr'	'MSE_tt' 'MSE_all'	'AIC_tr' 'AIC_tt'	'AIC_all' 'VAR' 'npar' 'time(s)'};
        table_dnn{2,1}=hnet.tag;table_dnn{count,2}=hnet.dropout; table_dnn{count,3}=hnet.minibatch;
        table_dnn{count,4}=num2str(hnet.alfa0,'%.1e'); table_dnn{count,5}=num2str(hnet.alfaf,'%.1e'); table_dnn{count,6}=hnet.nter/hnet.nw;
                table_dnn{count,7}=[num2str(mean(runs.RMSEtr),'%.2e') ' (±' num2str(std(runs.RMSEtr),'%.2e') ')'];
        table_dnn{count,8}=[num2str(mean(runs.RMSEtt),'%.2e') ' (±' num2str(std(runs.RMSEtt),'%.2e') ')'];
        table_dnn{count,9}=[num2str(mean(runs.RMSEal),'%.2e') ' (±' num2str(std(runs.RMSEal),'%.2e') ')'];
        table_dnn{count,10}=[num2str(mean(runs.MSEtr),'%.2e') ' (±' num2str(std(runs.MSEtr),'%.2e') ')'];
        table_dnn{count,11}=[num2str(mean(runs.MSEtt),'%.2e') ' (±' num2str(std(runs.MSEtt),'%.2e') ')'];
        table_dnn{count,12}=[num2str(mean(runs.MSEal),'%.2e') ' (±' num2str(std(runs.MSEal),'%.2e') ')'];
        table_dnn{count,13}=[num2str(mean(runs.AICctr),'%.3g') ' (±' num2str(std(runs.AICctr),'%.3g') ')'];
        table_dnn{count,14}=[num2str(mean(runs.AICctt),'%.3g') ' (±' num2str(std(runs.AICctt),'%.3g') ')'];
        table_dnn{count,15}=[num2str(mean(runs.AICcal),'%.3g') ' (±' num2str(std(runs.AICcal),'%.3g') ')'];
        temp=0;
        for ic=1:data.nbatch
            temp=temp+sum(sum(data.batch(ic).y.^2./(data.ystd(:,ic)'.^2)))/((size(data.batch(ic).y,1)-1)*size(data.batch(ic).y,2)*data.nbatch);
        end
        table_dnn{count,16}=[num2str(100-mean(runs.MSEal)/temp*100,'%.3g') ' (±' num2str(std(runs.MSEal)/temp,'%.3g') ')'];
        table_dnn{count,17}=[num2str(hnet.nw,'%.4g')];
        table_dnn{count,18}=[num2str(mean(runs.time(end,:)),'%.3g') ' (±' num2str(std(runs.time(end,:)),'%.3g') ')'];
        if i==1
            results=[results;table_dnn];
        else
            results=[results;table_dnn(2:end,:)];
        end
    end
    namefile='structures_fit_results'; % name to save the excel file
    writecell(results,[namefile '.xlsx'])
    fprintf('--------------------------------------!!!!--------------------------------------------\n');
    fprintf('The mat files of each hybrid structure were save as names given in variable structures{i,2} \n');
    fprintf('The overall statistical results were saved in the excel file : %s.xlsx \n',namefile);
    fprintf('--------------------------------------------------------------------------------------\n');
end


