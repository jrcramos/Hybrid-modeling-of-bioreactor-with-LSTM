%%% ReadBatchData.m reads data.xlsx and preprocess the conc, feeds, volume
%%% CalcReactedAmount.m calculated the reacted amount and the amounts in the
%   bioreactor( data.mr and data.accum)
%%% RatesEstimation estimated the reactions rates mM/gDW that can be used
%   for comparison with rates estimated using the hybrid model


clc
clear all
close all
w = warning ('off','all');

[batch1]=ReadBatchData();save batch1.mat batch1
load batch1.mat 
[batch2]=CalcReactedAmount(batch1); save batch2.mat batch2
load batch2.mat 
[batch3]=RatesEstimation(batch2,0,1:9); save batch3.mat batch3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% mAb production rates are in mg/gDCW  batch3.rates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load batch3.mat 
data_cho=struct();

%%% convert biomass do gDW
dw=[216.1 287.496]*1e-12*1e6; % dry weight g/Mcell
tdw2=batch3.age>=100;DW=ones(numel(batch3.age),1)*dw(1);DW(tdw2)=dw(2); % mass during exponteial growth vs stationary phase
batch3.m_r_vcd=batch3.m_r_vcd.*repmat(DW,1,batch3.nbatch);
batch3.vcd_in.val=batch3.vcd_in.val.*repmat(DW,1,batch3.nbatch); % dDW
batch3.vcd.val=batch3.vcd.val.*repmat(DW,1,batch3.nbatch)*1e3; % gDW/L

for i=1:9
    data_cho(i).batchid=batch3.batchid(i);
    data_cho(i).time=batch3.age';
    data_cho(i).m_r=[batch3.m_r_vcd(:,i) batch3.m_r_product(:,i)]; 
    temp=[];mets=[];conc=[batch3.vcd.val(:,i) batch3.product.val(:,i)];
    accum=[batch3.vcd_in.val(:,i) batch3.product_in.val(:,i)];
    for ib=1:23
        temp=[temp, batch3.met_m_r(ib).val(:,i)];mets=[mets batch3.met_m_r(ib).name];
        conc=[conc batch3.met(ib).conc(:,i)];
        accum=[accum batch3.met_in(ib).val(:,i)];
    end
    data_cho(i).m_r=[data_cho(i).m_r temp]';
    
    data_cho(i).vol=batch3.samples.reactorvolume(:,i)';
    data_cho(i).conc=conc';
    data_cho(i).accum=accum';
    data_cho(i).val=batch3.rates(i).val; 
    data_cho(i).std=batch3.rates(i).std; 
    data_cho(i).names=batch3.rates(i).names'; data_cho(i).names(1)={'Xv'};
end

data=data_cho;
save data.mat data