function [batch]=RatesEstimation(batch,iplot,varargin)
if nargin==3
    batchs=varargin{1};
else
    batchs=1:batch.nbatch;
end
fprintf('___________________________________________________________________________\n\n')
fprintf('Rates estimation %2.2f%% concluded.\n', 0)
% % find induction time;
% alli=[];
% for i=1:9
%      ids=batch.perbatch(i).ids;
%     all=batch.age_all;
%     age=batch.age;
%     x=batch.vcd.val(:,i);lac=batch.met(13).conc(:,i);
%     id=find(age==96);id=id(1);
%     diffs=all-age(id);diffs(diffs<0)=inf;
%     [~,id_in]=min(diffs);
%     diffs=ids-id;diffs(diffs>0)=-inf;
%     [~,id_ind]=max(diffs);
%     alli=[alli;id_ind-1];
%     figure('Name',['VCD_' num2str(i)]);
%     plot(age(~isnan(x)),x(~isnan(x)),'o');
%     lac=lac*(nanmax(x)/nanmax(lac));
%     hold on; plot(age(~isnan(lac)),lac(~isnan(lac)),'ro');
%     clc
%     a=[age(ids) (1:numel(ids))']'
% end
colors=cool(48);

p=[ 0,36, 38,84, 100,144, 146,240; %1
    0,36, 38,84, 100,144, 146,240; %2
    0,36, 38,84, 100,144, 146,240; %3
    0,36, 38,84, 100,144, 146,240; %4
    0,36, 38,84, 100,144, 146,240; %5
    0,36, 38,84, 100,144, 146,240; %6
    0,36, 38,84, 100,144, 146,240; %7
    0,36, 38,84, 100,144, 146,240; %8
    0,36, 38,84, 100,144, 146,240; %9
    ];

phases=p;
for i=1:9
    ids=batch.perbatch(i).ids;age=batch.age;
    A = repmat(age,[1 length(p(i,:))]);
    [~,temp] = min(abs(A-p(i,:)));
    A = repmat(ids,[1 length(temp)]);
    [~,temp] = min(abs(A-temp));
    phases(i,:)=temp;
end
ic=[1 12 24 48];
rates=struct();
dw=[216.1 216.1 287.496 287.496]*1e-12*1e6; % dry weight g/Mcell
names_met=[{'VCD'}; {'mAb'}];
for i= batchs
    close all
    rates(i).batchid=batch.batchid(i);
    
    ids=batch.perbatch(i).ids;
    
    ivc=NaN(numel(batch.age),1);ivc(ids)=batch.ivc(i).val;
    vcd_mr=batch.m_r_vcd(:,i);
    if iplot
        figure('Name','VCD');
        plot(ivc(~isnan(vcd_mr)),vcd_mr(~isnan(vcd_mr)),'o','MarkerSize',20,'MarkerEdgeColor',[0 0 0],'Color',[0 0 0],'HandleVisibility','off')
        hold on
        for ki=find((~isnan(vcd_mr)))
            text(ivc(ki), vcd_mr(ki)*1.075, string(batch.vcd.val(ki,i)),'FontSize',9);
            hold on
        end
        ylim([nanmin(vcd_mr)*1.1 nanmax(vcd_mr)*1.1]); xlim([-0.5e5 nanmax(ivc)*1.1])
    end
    age=batch.age;
    ivc_all=batch.ivc_all(i).val';
    vcd_mr=batch.m_r_vcd(:,i); vcd_mr=interp1(age(~isnan(vcd_mr)),vcd_mr(~isnan(vcd_mr)),batch.age_all,'linear','extrap');
    volume=batch.samples.reactorvolume(:,i)';volume=interp1(batch.age,volume,batch.age_all,'linear','extrap');
    std_vc=std(batch.vcd.val(:,i),'omitnan');
    for ip=1:4
        idphase=batch.age(ids(phases(i,[ip*2-1 ip*2-1+1])));
        A = repmat(batch.age_all',[1 length(idphase)]);[~,ids_] = min(abs(A-idphase'));

        vcd_mr_phase=vcd_mr(ids_(1):ids_(2)); ivc_phase=ivc_all(ids_(1):ids_(2));volume_phase=volume(ids_(1):ids_(2));
        vcd_mr_phase=vcd_mr_phase+std_vc*volume_phase.*randn(size(vcd_mr_phase))*0.1;% add error to generated points
        
        mod=fitlm(ivc_phase, vcd_mr_phase,'linear','RobustOpts','on');
        ci=coefCI(mod); ci = (ci(2,2)-ci(2,1))/2; val = mod.Coefficients.Estimate(2);
        rates(i).val(1,ip)=val;rates(i).std(1,ip)=ci;
        if iplot
            plot(ivc_phase(1),vcd_mr_phase(1),'-o','MarkerSize',7.5,'MarkerFaceColor',colors(ic(ip),:),'MarkerEdgeColor',colors(ic(ip),:),'HandleVisibility','on')
            plot(ivc_phase,vcd_mr_phase,'-o','MarkerSize',7.5,'MarkerFaceColor',colors(ic(ip),:),'MarkerEdgeColor',colors(ic(ip),:),'HandleVisibility','off')
        end
    end
    legend('P1','P2','P3','P4', 'location','best')

    if iplot
        xlabel('gDW x hour','fontsize',18);ylabel('Mcell','fontsize',18)
        title(['VCD ',batch.batchid{i}],'FontName','Times','fontsize',18)
    end
    
    
    mAb_mr=batch.m_r_product(:,i);
    if iplot
        figure('Name','mAb');
        plot(ivc(~isnan(mAb_mr)),mAb_mr(~isnan(mAb_mr)),'o','MarkerSize',20,'MarkerEdgeColor',[0 0 0],'Color',[0 0 0],'HandleVisibility','off')
        hold on
        ylim([nanmin(mAb_mr)*1.1 nanmax(mAb_mr)*1.1]); xlim([-0.5e5 nanmax(ivc)*1.1])
    end
    
    mAb_mr(age==0)=0;mAb_mr(age==24)=0;mAb_mr(age==60)=0;
    try 
        mAb_mr=interp1(age(~isnan(mAb_mr)),mAb_mr(~isnan(mAb_mr)),batch.age_all,'linear','extrap');
    catch kl
        1;
    end
    std_mAb=std(batch.product.val(:,i),'omitnan');
    for ip=1:4
        idphase=batch.age(ids(phases(i,[ip*2-1 ip*2-1+1])));
        A = repmat(batch.age_all',[1 length(idphase)]);[~,ids_] = min(abs(A-idphase'));

        mAb_mr_phase=mAb_mr(ids_(1):ids_(2)); ivc_phase=ivc_all(ids_(1):ids_(2));volume_phase=volume(ids_(1):ids_(2));
        mAb_mr_phase=mAb_mr_phase+std_mAb.*volume_phase.*randn(size(mAb_mr_phase))*1e-3*0.05;% add error to generated points
        try
            mod=fitlm(ivc_phase*dw(ip), mAb_mr_phase,'linear','RobustOpts','on');
        catch kl
            1;
        end
        ci=coefCI(mod); ci = (ci(2,2)-ci(2,1))/2; val = mod.Coefficients.Estimate(2);
        rates(i).val(2,ip)=val;rates(i).std(2,ip)=ci;
        if ip==1
            rates(i).val(2,ip)=0;
        end
        if iplot
            plot(ivc_phase,mAb_mr_phase,'-o','MarkerSize',10,'MarkerFaceColor',colors(ic(ip),:),'MarkerEdgeColor',colors(ic(ip),:),'HandleVisibility','off')
        end
    end
    
    if iplot
        xlabel('Mcell x hour','fontsize',18);ylabel('mM','fontsize',18)
        title(sprintf('mAb ',batch.batchid {i}),'FontName','Times','fontsize',18)
    end
    
    
    for k=1:size(batch.met,2) 
       met_mr=batch.met_m_r(k).val(:,i)'; 
       
       if iplot
           figure('Name',batch.met_m_r(k).name{:});
           plot(ivc(~isnan(met_mr)),met_mr(~isnan(met_mr)),'o','MarkerSize',20,'MarkerEdgeColor',[0 0 0],'Color',[0 0 0],'HandleVisibility','off')
           hold on
           ylim([nanmin(met_mr)*1.1 nanmax(met_mr)*1.1]); xlim([-0.5e5 nanmax(ivc)*1.1])
       end
       try
           met_mr_all=interp1(age(~isnan(met_mr)),met_mr(~isnan(met_mr)),batch.age_all,'linear','extrap');
       catch
           met_mr_all=nan(1,numel(batch.age_all));
       end
       std_met=std(batch.met_m_r(k).val(:,i),'omitnan');
       for ip=1:4
           flag=0;
           try
               idphase=batch.age(ids(phases(i,[ip*2-1 ip*2-1+1])));
               A = repmat(batch.age_all',[1 length(idphase)]);[~,ids_] = min(abs(A-idphase'));
               
               met_mr_phase=met_mr_all(ids_(1):ids_(2)); ivc_phase=ivc_all(ids_(1):ids_(2));volume_phase=volume(ids_(1):ids_(2));
%                met_mr_phase=met_mr_phase+met_mr_phase.*volume_phase./max(volume).*randn(size(met_mr_phase))*0.075;% add error to generated points
               met_mr_phase=met_mr_phase+std_met.*volume_phase.*randn(size(met_mr_phase))*1e-3*0.2;% add error to generated points
               
  
               mod=fitlm(ivc_phase*dw(ip), met_mr_phase,'linear','RobustOpts','on');
               ci=coefCI(mod); ci = (ci(2,2)-ci(2,1))/2; val = mod.Coefficients.Estimate(2);
               rates(i).val(2+k,ip)=val;rates(i).std(2+k,ip)=ci;flag=1;
           catch kl
               rates(i).val(2+k,ip)=NaN;rates(i).std(2+k,ip)=NaN;
           end
           if iplot && flag==1
               plot(ivc_phase,met_mr_phase,'-o','MarkerSize',10,'MarkerFaceColor',colors(ic(ip),:),'MarkerEdgeColor',colors(ic(ip),:),'HandleVisibility','off')
           end
       end
        if iplot && flag==1
            xlabel('Mcell x hour','fontsize',18);ylabel('mmol','fontsize',18)
            title(sprintf(batch.met_m_r(k).name{:}, ' ',batch.batchid {i}),'FontName','Times','fontsize',18)
        end
        try names_met(k+2)=batch.met_m_r(k).name; catch  names_met(k+2)={batch.met_m_r(k).name};end
    end
    disp([repmat(char(8), 1, 36)])
    fprintf('Rates estimation %2.2f%% concluded.\n', i/9*100)
    
end


for i=1:9
    rates(i).names=names_met;
end


batch.rates=rates;

end