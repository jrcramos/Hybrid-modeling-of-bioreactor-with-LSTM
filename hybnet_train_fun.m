function []=hybnet_train_fun(simfile,layers,Indtr,Indcr,vargin)
rng default

ymode=2; % 1 : y is reacted mass ; 2 : y is concentration
dopca=1; % 1: do PCA of reacted amount; 2 : identity matrix is generated instead of doing PCA
npcs=4; % number of principal compoenents to uso on PCA

if nargin==1
    load([simfile '.mat'])
    dotraining=0; % 1: run model fit <-----> 0 : run simulation of simfile given
    niter=[];nruns=[];npcs=4;
else
    niter=vargin{1}; %number of iteration for ADAM
    nruns=vargin{2}; % number of training repetitions, up to 200 pre set permutattions
    dotraining=1; % 1: run model fit <-----> 0 : run simulation of given
    npcs=vargin{3};% number of principal compoenents to uso on PCA
    % Create and test a hybrid net............  
    hnet = hnetcreate(layers);
    s='In(27)';
    for il=1:hnet.nl
        if strcmp(hnet.layers{1, il}.type ,'LSTM')
            s=[s '-LSTM' '(' num2str(hnet.layers{1, il}.ny) ')'];
        else
            temp=hnet.layers{1, il}.node;temp(1)=upper(temp(1));
            s=[s '-' temp '(' num2str(hnet.layers{1, il}.ny) ')'];
        end
    end
    s=[s '-Rate'  '(' num2str(hnet.layers{1, il}.ny) ')' '-ODE(25)'];
    s=replace(s,'Linear','Lin');s=replace(s,'Relu','ReLU'); % structure name
    hnet = hnetcreate(layers,s);
end
if dotraining
    fprintf('-------------------------!!!!-----------------------------------\n');
    fprintf('Option do model training activated \n');
    fprintf('----------------------------------------------------------------\n');
else
    fprintf('--------------------------------------!!!!--------------------------------------------\n');
    fprintf('Option do model training not activated, doing model simulation using simfile: %s \n',simfile);
    fprintf('---------------------------------------------------------------------------------------\n');
end
    
[data]=readdata(ymode,dopca,npcs);






%do training
runs=struct('indtr',[],'wfinal',[],'RMSEtr',[],'RMSEtt',[],'RMSEal',[],...
    'MSEtr',[],'MSEtt',[],'MSEal',[],'AICctr',[],'AICctt',[],'AICcal',[],'itern',[],...
    'time',[],'mse',[],'mse2',[],'hnet',{},'tt',{},'data_iter',{},'tt2',{},'mse3',{},'mse4',{},...
    'crossv',[],'MSEcr',[],'RMSEcr',[]);


hnet.ymode=ymode;
hnet.dopca=dopca;
hnet.minibatch = 1; % minibatch size ]0-1], 1: no minibatch is used
hnet.dropout =  0.1; % regularization stochastic weights dropout [0-1[, 0: no weights are randomly set to zero
hnet.alfa=1e-3;  %ADAM learining rate
hnet.alfa0=hnet.alfa; % inital alfa in linear decay
hnet.alfaf=hnet.alfa;%0.0005; % final alpha in linear decay

hnet.nter=niter;
simfile=simfile;% name to save results


resume=0;
rng(1803199136)
hnet=hnetinitw(hnet);

runsv=1:nruns;
if resume
   load runs.mat
    if runs.irun==nruns && dotraining
        fprintf('increase the Nruns \n')
        return
    else
        runsv=runs.irun+1:nruns; hnet=runs.hnet{end};
    end
    
end


if dotraining
    %load('Indtr.mat');load('Indcr.mat');
    if resume
        runs.indtr(:,runsv)=Indtr(runsv,:)';
        if ~isempty(Indcr)
            runs.indcr(:,runsv)=Indcr(runsv,:)';
        end
    else
        runs(1).indtr(:,1:nruns)=Indtr(1:nruns,:)';
        if ~isempty(Indcr)
            runs(1).indcr(:,1:nruns)=Indcr(1:nruns,:)';
        end
    end
    rng(1803199136+1)
    hnet=hnetinitw(hnet);[w0]=hnetgetw(hnet);
    for irun=runsv        
        % update training and testing batchs
        for i=1:data.nbatch
            data.batch(i).istrain=2;
        end
        for i=1:length(Indtr(1,:))
            data.batch(Indtr(irun,i)).istrain=1;
        end
        if ~isempty(Indcr)
            for i=1:length(Indcr(1,:))
                data.batch(Indcr(irun,i)).istrain=3;
            end
        end
        
        hnet=hnetsetw(hnet,w0);
        [hnet,mse, iters, time,tt,mse2,tt2,mse3,mse4,crossv]=hnettrain_new_v1(hnet,data,niter,1+nruns-irun);
%         [hnet,mse, iters, time]=arrayfun(@hnettrain,hnet,data,niter,1+nruns-irun);
        
        %final simulation
        [hnet,data,mse_al,mse_tr,mse_tt,rmse_al,rmse_tr,rmse_tt,AICc_al,AICc_tr,AICc_tt,mse_cr,rmse_cr]=hnetsimul(hnet,data,0,1);
        [wfinal]=hnetgetw(hnet);
        
        
        runs.wfinal(:,end+1)=wfinal;
        runs.RMSEtr(:,end+1)=rmse_tr;
        runs.RMSEtt(:,end+1)=rmse_tt;
        runs.RMSEal(:,end+1)=rmse_al;
        runs.MSEtr(:,end+1)=mse_tr;
        runs.MSEtt(:,end+1)=mse_tt;
        runs.MSEal(:,end+1)=mse_al;
        runs.AICctr(:,end+1)=AICc_tr;
        runs.AICctt(:,end+1)=AICc_tt;
        runs.AICcal(:,end+1)=AICc_al;
        runs.itern(:,end+1)=iters;
        runs.time(:,end+1)=time;
        runs.mse(:,end+1)=mse;
        runs.mse2(:,end+1)=mse2;
        runs.hnet(end+1)={hnet};
        runs.data=data;
        runs.irun=irun;
        runs.tt(end+1)={tt};
        runs.tt2(end+1)={tt2};
        runs.mse3(end+1)={mse3};
        runs.mse4(end+1)={mse4};
        runs.data_iter(end+1)={data};
        runs.crossv(:,end+1)=crossv;
        runs.MSEcr(:,end+1)=mse_cr;
        runs.RMSEcr(:,end+1)=rmse_cr;
        
        fprintf('--------------------------\\--------------------------------\n')
        fprintf('Finished run %d out of %d  \n', irun, nruns)
        fprintf('---------------------------\\-------------------------------\n\n')
        save runs.mat runs
    end
    
    
    % find best run and prepare for plots
    [ind]=get_best_network(runs);
    wfinal=runs.wfinal(:,ind);
    for i=1:data.nbatch
        data.batch(i).istrain=2;
    end
    Ind=runs.indtr(:,ind);
    for i=1:numel(Ind)
        data.batch(Ind(i)).istrain=1;
    end
    
    if ~isempty(Indcr)
        Indcr=runs.indcr(:,ind);
        for i=1:numel(Indcr)
            data.batch(Indcr(i)).istrain=3;
        end
    end
    % % Simulation and plot ###
    iplot=1;
    hnet=runs.hnet{ind};hnet.bestperm=ind;
    
    fprintf('--------------------------\\--------------------------------\n')
    fprintf('Best permuation was %d out of %d  \n', ind, nruns)
    fprintf('---------------------------\\-------------------------------\n')
    [~]=hnetsimul(hnet,data,iplot,1);
    
    save([simfile '.mat'], 'hnet', 'data','wfinal','runs','nruns','ind')
    
    
else
    load([simfile '.mat']);
    iplot=1;
    fprintf('--------------------------\\--------------------------------\n')
    fprintf('Best permuation was %d out of %d  \n', ind, nruns)
    fprintf('---------------------------\\-------------------------------\n')
    [~]=hnetsimul(hnet,data,iplot,1);
    [hnet,data,mse_al,mse_tr,mse_tt,rmse_al,rmse_tr,rmse_tt,AICc_al,AICc_tr,AICc_tt]=hnetsimul(hnet,data,0,1);
end


legends={};
for ib=1:nruns;legends(end+1)={['run#' num2str(ib)]};end
allids=1:runs.irun;
if runs.irun>1; allids(ind)=[]; end

figure('Name','MSE','units','normalized','outerposition',[0 0.1 0.5 0.4])
% subplot(2,1,1)
% semilogy(runs.itern(:,allids),(runs.mse(:,allids)),'Color','#00BFFF','HandleVisibility','off');hold on;
% semilogy(runs.itern(:,allids(1)),(runs.mse(:,allids(1))),'Color','#00BFFF','LineWidth',1.5,'HandleVisibility','on')
% semilogy(runs.itern(:,ind),(runs.mse(:,ind)),'Color','#000000','LineWidth',1.5,'HandleVisibility','on');
% ax = gca;ax.TickLength = [0.02, 0.02];ax.LineWidth = 100*0.01;
% legend({'other permutation' 'best permutation'},'Location','best')
% xlim([0 max(max(runs.itern))*1.02]);ylim([0 max(max(runs.mse))*1.2]);
% % legend(legends,'Location','bestoutside')
% xlabel('#iteration','FontSize',12,'FontWeight','Bold','FontName','Arial')
% ylabel('MSE','FontSize',12,'FontWeight','Bold','FontName','Arial')
% % text(runs.itern(end,1)*0.95,runs(1).mse(1,1)*0.6*1.1,'A','FontSize',12,'FontWeight','Bold','FontName','Arial')
% set(gca,'FontSize',12,'FontName','Arial')
% subplot(2,1,2)
semilogy(runs.itern(:,allids),(runs.mse(:,allids)),'Color','#00BFFF','HandleVisibility','off');hold on;
semilogy(runs.itern(:,allids(1)),(runs.mse(:,allids(1))),'Color','#00BFFF','LineWidth',1.5,'HandleVisibility','on');hold on;
semilogy(runs.itern(:,ind),(runs.mse(:,ind)),'Color','#000000','LineWidth',1.5,'HandleVisibility','on')
semilogy(runs.itern(:,allids),(runs.crossv(:,allids)),'Color','r','HandleVisibility','off');hold on;
semilogy(runs.itern(:,allids(1)),(runs.crossv(:,allids(1))),'Color','r','LineWidth',1.5,'HandleVisibility','on');hold on;
semilogy(runs.itern(:,ind),(runs.crossv(:,ind)),'Color','b','LineWidth',1.5,'HandleVisibility','on')

ax = gca;ax.TickLength = [0.02, 0.02];ax.LineWidth = 100*0.01;
xlim([0 max(max(runs.itern))*1.02]);ylim([0 max(max(runs.mse))*1.2]);
legend({'MSE train' ' MSE train best' 'MSE crossv' ' MSE crossv best'},'Location','best')
xlabel('#iter','FontSize',12,'FontWeight','Bold','FontName','Arial')
ylabel('MSE','FontSize',12,'FontWeight','Bold','FontName','Arial')
% text(runs.time(end,1)*0.975,runs(1).mse(1,1)*0.6*1.1,'B','FontSize',12,'FontWeight','Bold','FontName','Arial')
set(gca,'FontSize',12,'FontName','Arial')


plot_rates(hnet,data)
plot_train_test(hnet,data)


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ind]=get_best_network(runs)
rmseall=sum([runs.RMSEtr ;runs.RMSEtt]);
ind = find(rmseall==min(rmseall));ind=ind(1);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=plot_rates(model,data)
figure('Name','Rates','units','normalized','outerposition',[0.07 0 0.4 0.55])
data.ylabel(1) = '\mu';
% tiledlayout(7,4,'Padding','none','TileSpacing','compact')

for k=1:numel(data.ylabel)
    subplot(6,5,k)  
	hold on
    for ib=1:data.nbatch
		x = data.age_h  ;
		y = data.batch(ib).rates(:,k);
		if data.batch(ib).istrain==1
			trainingPlot = plot(x,y,...
				'MarkerSize',5,'Color','#0080ff','LineStyle','-');
		else
			testPlot = plot(x,y,...
				'MarkerSize',5,'Color','#2bbf5c','LineStyle','-');
        end
    end
    
    if k>1;ylabel('mmol/gDW/h');else;ylabel('h^{-1}');end
    xlabel('time (h)');
    set(gca,'FontSize',6);
    set(gca,'FontWeight','bold');
    set(gca,'LineWidth',2);
    title(data.ylabel(k))
    axis square
    grid on
end
lgd = legend([trainingPlot testPlot],{'Training', 'Testing'});%,'Location','bestoutside');
lgd.FontName = "Liberation Sans";
lgd.FontSize = 10;
lgd.FontWeight = 'bold';
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=plot_train_test(model,data)
figure('Name','fit individually','units','normalized','outerposition',[0.15 0.12 0.4 0.55])
% tiledlayout(7,4,'Padding','none','TileSpacing','compact')


for k=1:numel(data.ylabel)
    subplot(6,5,k)  
	hold on
    for ib=1:data.nbatch
		x = data.batch(ib).y(:,k);
		y = data.batch(ib).ysim(:,k);
        if data.batch(ib).istrain==1
            trainingPlot = plot(x,y,...
                'MarkerSize',5,'Color','#0080ff','Marker','.','LineStyle','none');
        else
            testPlot = plot(x,y,...
                'MarkerSize',5,'Color','#2bbf5c','Marker','x','LineStyle','none');
        end
    end
    
    xlabel('experimental');ylabel('prediction');
    set(gca,'FontSize',6);
    set(gca,'FontWeight','bold');
    set(gca,'LineWidth',2);
    title(data.ylabel(k))
    axis square
    grid on

end

lgd = legend([trainingPlot testPlot],{'Training', 'Testing'});%,'Location','bestoutside');
lgd.FontName = "Liberation Sans";
lgd.FontSize = 10;
lgd.FontWeight = 'bold';

end

function [data]=readdata(ymode,dopca,npcs)
% load data
[x,y,accum,formed,vol,age_h,xlabels,ylabels,reacted_masses]=load_data_hybmod(ymode);

xIdx = 1;   % biomass column gDW/L

xlabels(end+1:end+2)={'Vol' 'time'};


% Data preprocessing ####
deltat=age_h(2)-age_h(1);
xall =[];
yall =[];
ymaxall=[];
for ib=1:length(x)
    xbatch = cell2mat(x(ib));v=vol{ib};
    xb=[xbatch; v; age_h ];%volume and time may considered in the input
    data.flag_t_v=1; % if volume and time are considered in the x
    xall=[xall, xb];
    ybatch = cell2mat(y(ib));
    yall=[yall, ybatch(:,:)];ymaxall(:,ib)=max(ybatch(:,:)');
    data.batch(ib).x=xb';
    data.batch(ib).y=ybatch(:,:)';
    data.batch(ib).accum=cell2mat(accum(ib));
    data.batch(ib).formed=cell2mat(formed(ib));
    data.batch(ib).vol=cell2mat(vol(ib));
    data.batch(ib).mx=xbatch(xIdx,:).*v/1000*deltat;
    data.batch(ib).istrain=2; %set for testing
    %   calaculate the mean cell mass in the time interval
    % 	for i=1:length(data.batch(ib).mx)-1
    % 		if data.batch(ib).mx(i)~=data.batch(ib).mx(i+1)
    % 			data.batch(ib).mx(i)=(data.batch(ib).mx(i+1)-data.batch(ib).mx(i))/...
    % 				log(data.batch(ib).mx(i+1)/data.batch(ib).mx(i));
    % 		end
    % 	end
end

% do PCA of reacted mass
if dopca==1
    reacted_masses_max=max(abs(reacted_masses));
    reacted_masses_pca=reacted_masses;
    for i=1:size(reacted_masses,1)
        reacted_masses_pca(i,:)=reacted_masses(i,:)./reacted_masses_max;
    end
    nscores=npcs;
    [coeff,scores,~,~,explained] = pca(reacted_masses_pca,'NumComponents',nscores,'Algorithm','svd','Centered',false);
    data.coeff=coeff;
    data.scaling=reacted_masses_max;vari=cumsum(explained);
    fprintf('PCA of reacted amounts performed cuptured %.2f variance of data \n',vari(nscores));
else
    data.coeff=diag(ones(1,size(yall,1)));
    data.scaling=ones(1,size(yall,1));
end

xmean = mean(xall,2);
ymax = max(abs(yall)')';
xstd = std(xall,[],2);
ymean = mean(yall,2);
ystd =std(yall,[],2);
data.xmean=xmean';
data.xmin=min(xall,[],2)';
data.xmax=max(xall,[],2)';
data.xstd=xstd';
data.ymean=ymean;
data.ystd=repmat(ymax,1,20);%repmat(ystd,1,20);%ymaxall;%ystd;
data.ystd2=data.ystd;
% data.ystd=data.ystd.*weights; %higher importance some conc such as, Glc, Lac, pyr
data.cnoise=[0.1 ones(1,26)*0.2]'; % vector of noise to add model input (concentrations), when fitting reacted masses, ymode=1, 0.1 -> 10%
data.ymax=ymax;
data.nbatch=length(x);
data.deltat=deltat;
data.ylabel=ylabels;
data.xlabel=xlabels;
data.age_h=age_h;

% create data with noise
% rng(1803199136)
% noise=data.cnoise;
% for ib=1:data.nbatch
%     for it=1:size(data.batch(ib).y,1)
%         data.batch(ib).ynoise(it,:)= data.batch(ib).y(it,:).*(1+randn(numel(noise),1).*noise)';
%     end
% 
% end


xnorm = xall';
max_=max(xall');
for i=1:size(xnorm,2)
    if 1
        data.xnorm=1;
        xnorm(:,i)=(xnorm(:,i)-xmean(i))/xstd(i);
    elseif 0
        data.xnorm=2;
        xnorm(:,i)=(xnorm(:,i)-min(xnorm(:,i)))./(max(xnorm(:,i))-min(xnorm(:,i)));    
    else
        data.xnorm=3;
    end
end

age_end=size(xb,2)-1;
count=1;
for ib=1:data.nbatch
    data.batch(ib).x = xnorm(count:count+age_end,:);
    
    count=count+age_end+1;
end


% save data_temp.mat data

end