function [hnet,data,mse_al,mse_tr,mse_tt,rmse_al,rmse_tr,rmse_tt,AICc_al,AICc_tr,AICc_tt,mse_cr,rmse_cr]=hnetsimul(hnet,data,iplot,normres) %##########################
if nargin<3
    iplot=1;
    normres=0;
elseif nargin<4
    normres=0;
end
mse_tr = 0;
count_tr = 0;
mse_cr = 0;
count_cr = 0;
mse_tt = 0;
count_tt = 0;
mse_al = 0;
count_al = 0;

if iplot==1
    figure('Name','Overall Fit','units','normalized','outerposition',[0.3 0.1 0.5*0.5 0.7*0.6])
end
%#########################################################################
% ?????????????????????????????????????????????????????????????????????????
% activate only to predicited concentrations based on model fitted to reacted mass
%     trainid=[data.batch.istrain];
%     load data_temp
%     hnet.ymode=2;
%     for i=1:data.nbatch
%         data.batch(i).istrain=trainid(i);
%     end
% ?????????????????????????????????????????????????????????????????????????
%###########################################################################


minmaxx=[0 0];minmaxy=[0 0];
allpred=[];allsim=[];allcross=[];
try hnet.ymode; catch; hnet.ymode=1; end%Joao Ramos
stoichiometry=data.coeff.*data.scaling';
if data.flag_t_v
    data.xmean=data.xmean(1:numel(data.ylabel));data.xstd=data.xstd(1:numel(data.ylabel));
    data.xmax=data.xmax(1:numel(data.ylabel));data.xmin=data.xmin(1:numel(data.ylabel));
end
for ib=1:data.nbatch
    hnet=hnet_init_state(hnet); 
    xbatch= data.batch(ib).x;
    ybatch= data.batch(ib).y;    
    
    mx= data.batch(ib).mx; %
    mr =  zeros(size(stoichiometry,1),1);%zeros(hnet.ny,1);
    if hnet.ymode~=1  %João Ramos
        predc=data.batch(ib).accum(1,:)/(data.batch(ib).vol(1)*1e-3);
    end

    data.batch(ib).ysim = zeros(size(ybatch));
     if hnet.ymode~=1  %João Ramos
         data.batch(ib).ysim(1,:)=predc;
     end
    data.batch(ib).rates = zeros(size(ybatch));
    if strcmp(hnet.layers{hnet.nl}.type,'Metabolic')
        data.batch(ib).rates_full = zeros(size(ybatch,1),numel(hnet.layers{hnet.nl}.v));
    end
    
    if data.flag_t_v
        time_volume=data.batch(ib).x(:,end-1:end); % t and V are fixed inputs
    end
    for t=2:size(ybatch,1)
        if data.flag_t_v
            t_v=time_volume(t-1,:);
        else
            t_v=[];
        end
       if hnet.ymode ==1 %João Ramos     
          xinp  = [xbatch(t-1,:)' t_v];
       else
           if data.xnorm==1
               xinp  =[[predc(t-1,:)]-data.xmean./(data.xstd) t_v];%João Ramos
           elseif data.xnorm==2
               % xinp  =(([predc(t-1,:) data.age_h(t-1)]-data.xmin)./(data.xmax-data.xmin))*2-1;%João Ramos
               xinp  =[(([predc(t-1,:) ]-data.xmin)./(data.xmax-data.xmin))*2-1 t_v];%João Ramos
           else
               fprintf('Input not normalized \n')
           end
       end
       hnet=hnetff(hnet,xinp);
       if hnet.ymode~=1  %João Ramos
           mx(t-1)= predc(t-1,1)*data.deltat*data.batch(ib).vol(t-1)/1000;
       end
       mr = mr + mx(t-1)*(stoichiometry*hnet.y);
       
       if hnet.ymode~=1  %João Ramos
           predc(t,:)=(data.batch(ib).accum(t,:)+mr')/(data.batch(ib).vol(t)*1e-3);
       end
       
       if hnet.ymode==1  %João Ramos
           data.batch(ib).ysim(t,:)=mr';
       else %João Ramos
           data.batch(ib).ysim(t,:)= predc(t,:);
       end
       data.batch(ib).rates(t,:) =stoichiometry*hnet.y;
       
       if strcmp(hnet.layers{hnet.nl}.type,'Metabolic')
           data.batch(ib).rates_full(t,:)=hnet.layers{hnet.nl}.v;
       end
%        v=data.batch(ib).rates_full(:,hnet.layers{hnet.nl}.Ir);
       
       ytarg = ybatch(t,:)';
       if hnet.ymode==1  %João Ramos
           if normres
               res = (ytarg-mr)./data.ystd2(:,ib); 
           else
               res = (ytarg-mr)./data.ystd(:,ib); 
           end
       else
           if normres
               res = (ytarg- predc(t,:)')./data.ystd2(:,ib); %João Ramos
           else
               res = (ytarg- predc(t,:)')./data.ystd(:,ib); %João Ramos
           end
       end
       
       count_al = count_al+numel(res);    
       mse_al = mse_al + res'*res;
       
       if data.batch(ib).istrain == 1%-------           
            count_tr = count_tr+numel(res);    
            mse_tr = mse_tr + res'*res;
       elseif data.batch(ib).istrain == 3%-------           
            count_cr = count_cr+numel(res);    
            mse_cr = mse_cr + res'*res;
       else
            count_tt = count_tt+numel(res);    
            mse_tt = mse_tt + res'*res;
       end
       
    end
    
    data.batch(ib).rates(1,:)=data.batch(ib).rates(2,:);
    if strcmp(hnet.layers{hnet.nl}.type,'Metabolic')
        data.batch(ib).rates_full(1,:)=data.batch(ib).rates_full(2,:);
    end
    
    
end
if iplot == 1
    

    for ib=1:data.nbatch %plot training
        normy=repmat(data.ymax',size(data.batch(ib).y,1),1);
        
        minmaxx(1)=min([minmaxx(1) min(data.batch(ib).y./normy)]);
        minmaxx(2)=max([minmaxx(2) max(data.batch(ib).y./normy)]);
        minmaxy(1)=min([minmaxy(1) min(data.batch(ib).ysim./normy)]);
        minmaxy(2)=max([minmaxy(2) max(data.batch(ib).ysim./normy)]);
        if data.batch(ib).istrain == 1
            allsim=[allsim;[data.batch(ib).y./normy data.batch(ib).ysim./normy]];
            trainingPlot = plot(data.batch(ib).y./normy,data.batch(ib).ysim./normy,...
                'Marker','.','LineStyle','none','MarkerSize',8,...
                'Color','#0080ff', 'HandleVisibility','off');
            hold on
            
        end
    end
    for ib=1:data.nbatch %plot test on top
        
        if data.batch(ib).istrain == 3
            allcross=[allcross;[data.batch(ib).y./normy data.batch(ib).ysim./normy]];
            crossvPlot = plot(data.batch(ib).y./normy,data.batch(ib).ysim./normy,...
                'Marker','.','LineStyle','none','MarkerSize',5.5,...
                 'Color','#44a6c6','HandleVisibility','off');
            hold on
        end
        
        if data.batch(ib).istrain == 2
            allpred=[allpred;[data.batch(ib).y./normy data.batch(ib).ysim./normy]];
            testPlot = plot(data.batch(ib).y./normy,data.batch(ib).ysim./normy,...
                'Marker','.','LineStyle','none','MarkerSize',5.5,...
                 'Color','#2bbf5c','HandleVisibility','off');
            hold on
        end
    end
end


if iplot
     allsimx=allsim(:,1:size(allsim,2)*0.5); allsimy=allsim(:,1+size(allsim,2)*0.5:end); 
     allsimx=reshape(allsimx,numel(allsimx),1);allsimy=reshape(allsimy,numel(allsimy),1);
     allpredx=allpred(:,1:size(allpred,2)*0.5); allpredy=allpred(:,1+size(allpred,2)*0.5:end); 
     allpredx=reshape(allpredx,numel(allpredx),1);allpredy=reshape(allpredy,numel(allpredy),1);
     allcrossx=allcross(:,1:size(allcross,2)*0.5); allcrossy=allcross(:,1+size(allcross,2)*0.5:end); 
     allcrossx=reshape(allcrossx,numel(allcrossx),1);allcrossy=reshape(allcrossy,numel(allcrossy),1);
     mod=fitlm(allsimx,allsimy,'linear','RobustOpts','off');
     R2_tr=mod.Rsquared;
     mod=fitlm(allpredx,allpredy,'linear','RobustOpts','off');
     R2_tt=mod.Rsquared;
      if ~isempty(allcross)
          mod=fitlm(allcrossx,allcrossy,'linear','RobustOpts','off');
          R2_cr=mod.Rsquared;
      end
      
     plot(allsim(1,1),allsim(1,2),'Marker','.','LineStyle','none','MarkerSize',14,'Color','#0080ff', 'HandleVisibility','on');
     plot(allpred(1,1),allpred(1,2),'Marker','.','LineStyle','none','MarkerSize',14,'Color','#2bbf5c', 'HandleVisibility','on');
     if ~isempty(allcross)
         plot(allpred(1,1),allcross(1,2),'Marker','.','LineStyle','none','MarkerSize',14,'Color','#44a6c6', 'HandleVisibility','on');
     end
     if ~isempty(allcross)
     legend(['training (R^{2}= ' num2str(R2_tr.Ordinary,'%.3f') ')'],['crossv (R^{2}= ' num2str(R2_cr.Ordinary,'%.3f') ')'],...
         ['test (R^{2}= ' num2str(R2_tt.Ordinary,'%.3f') ')'],'Location','southeast');legend boxoff
     else
         legend(['training (R^{2}= ' num2str(R2_tr.Ordinary,'%.3f') ')'],['test (R^{2}= ' num2str(R2_tt.Ordinary,'%.3f') ')'],'Location','southeast');legend boxoff
     end
    xlabel('Experimental','FontWeight','Bold');ylabel('Prediction','FontWeight','Bold')
    set(gca,'FontSize',12,'FontName','Arial');yticks([0:0.5:1]);yticklabels([0:0.5:1]);
    ylim([ min([minmaxx(1) minmaxy(1)])*1.2 max([minmaxx(2) minmaxy(2)])*1.2])
    xlim([ min([minmaxx(1) minmaxy(1)])*1.2 max([minmaxx(2) minmaxy(2)])*1.2])
    try patchline([ min([minmaxx(1) minmaxy(1)])*1.2 max([minmaxx(2) minmaxy(2)])*1.2], [ min([minmaxx(1) minmaxy(1)])*1.2 max([minmaxx(2) minmaxy(2)])*1.2],'edgecolor','k','linewidth',1,'edgealpha',0.1,'HandleVisibility','off');catch; end
    xlim([-0.35 1.35]);ylim([-0.35 1.35])
    line([ -0.35 1.35],[ -0.35 1.35],'color',[.35 0.35 0.35],'HandleVisibility','off');hold on;
    line([ -0.35 1.15],[ -0.15 1.35],'LineStyle','--','color',[.68 0.68 0.68],'HandleVisibility','off');hold on;
    line([ -0.15 1.35],[ -0.35 1.15],'LineStyle','--','color',[.68 0.68 0.68],'HandleVisibility','off');hold on;
end
mse_tr = mse_tr/max(count_tr,1);
mse_cr = mse_cr/max(count_cr,1);
mse_tt = mse_tt/max(count_tt,1);
mse_al = mse_al/max(count_al,1);

nw=hnet.nw;
rmse_cr = sqrt(mse_cr);
rmse_tr = sqrt(mse_tr);
AICc_tr = count_tr*log(mse_tr)+2*nw + 2*nw*(nw+1)/(count_tr-nw-1);
rmse_tt = sqrt(mse_tt);
AICc_tt = count_tt*log(mse_tt)+2*nw + 2*nw*(nw+1)/(count_tt-nw-1);
rmse_al = sqrt(mse_al);
AICc_al = count_al*log(mse_al)+2*nw + 2*nw*(nw+1)/(count_al-nw-1);

if iplot==1
    fprintf('\n\nRESULTS:');
    fprintf('RMSE training: %f\n',rmse_tr);
    fprintf('RMSE crossvalidation: %f\n',rmse_cr);
    fprintf('RMSE testing: %f\n',rmse_tt);
%     fprintf('RMSE all: %f\n\n',rmse_al);
    fprintf('MSE training: %f\n',mse_tr);
    fprintf('MSE crossvalidation: %f\n',mse_cr);
    fprintf('MSE testing: %f\n',mse_tt);
%     fprintf('MSE all: %f\n\n',mse_al);
    fprintf('AICc: %f\n',AICc_tr);
%     fprintf('AICc testing: %f\n',AICc_tt);
%     fprintf('AICc all: %f\n\n',AICc_al);
end
end

