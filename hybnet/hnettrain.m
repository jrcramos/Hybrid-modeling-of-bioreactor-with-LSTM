function [hnetbest,mse, iters, time,tt,mse2,tt2]=hnettrain(hnet,data,niter,nruns)
% data = structfun(@gpuArray, data, 'UniformOutput', false);
% hnet = structfun(@gpuArray, hnet, 'UniformOutput', false);
% 
% verify gradients
% hnet_verify_grads(hnet,data);
% pause
iters=(1:niter)'; mse=zeros(niter,1);time=zeros(niter,1);mse2=zeros(niter,1);

hnetbest=hnet;fbest=inf;bestpenalty=inf; %Joao Ramos
hnet.iter=0;
cpu0 = cputime; iterper=round(0.05*niter);ncross=round(0.01*niter);
alfa0=hnet.alfa0;tt=[];tt2=[];
for i=1:niter
    [hnet,loss,count,penalty,loss2]=lossfun(hnet,data);
    
    cpui = cputime-cpu0;
    predcpu=cpui/i*niter;remcpu=predcpu*(1-i/niter);
    loss = loss/max(count,1);loss2=loss2/max(count,1);
    penalty=penalty/max(count,1);not_changed=1;
    if 1 %loss<fbest %&& penalty<(1+rand)*bestpenalty
        fbest=loss;hnetbest=hnet;%Joao Ramos
%         [~,~,~,~,mse_tt,~,~,~,~,~,~]=hnetsimul(hnet,data,0,1);
%         tt=[tt;[i mse_tt]];not_changed=0;
%         bestpenalty=penalty; %Joao Ramos
    else
        hnet.alfa0=alfa0;
    end
    
    if rem(i,ncross)==0 || i==1 || i==niter
         [~,~,~,~,mse_tt,~,~,~,~,~,~]=hnetsimul(hnet,data,0,1);
          [~,~,~,~,mse_tt2,~,~,~,~,~,~]=hnetsimul(hnet,data,0,0);
          tt=[tt;[i mse_tt]];tt2=[tt2;[i mse_tt2]];
    end
    if rem(i,iterper)==0 || i==1 || i==niter  
        fprintf('Iter:%5u    Fval: %2.2e    Fbest: %2.2e     alpha:%1.2e     Time: %5.1f    ETA: %5.1f    Total: %5.1f   ETA %d runs: %3.2fh \n',i, loss,fbest,hnet.alfa, cpui,remcpu, predcpu,nruns,(predcpu*nruns)/3600);
    end
    mse(i)=loss;mse2(i)=loss2;time(i)=cpui;
    hnet=hnet_update_weights(hnet);
end
end



function [hnet,mse,count,penalty,mse2]=lossfun(hnet,data) %#############################

mse = 0;penalty=0;mse2=0;
count = 0;
hnet = hnet_reset_bp(hnet);

if strcmp(hnet.layers{hnet.nl}.type,'Metabolic')
    ir=hnet.layers{hnet.nl}.Ir;
    nsiirt = hnet.layers{hnet.nl}.nullsi(ir,:)';
end
stoichiometry=data.coeff.*data.scaling';

for ib=1:data.nbatch
    
    if data.batch(ib).istrain == 1%-------
        
        hnet=hnet_init_state(hnet);
        
        xbatch= data.batch(ib).x;
        ybatch= data.batch(ib).y;
        
        mx= data.batch(ib).mx;
        mr = zeros(size(stoichiometry,1),1);%zeros(hnet.ny,1);
        if hnet.ymode~=1  %João Ramos
            predc=data.batch(ib).accum(1,:)/(data.batch(ib).vol(1)*1e-3);
        end
        dmrdy = zeros(size(stoichiometry,1),1);%zeros(hnet.ny,1);
        
        for t=2:size(ybatch,1)
            if hnet.ymode==1  %João Ramos
                xinp  = xbatch(t-1,:)';
                xinp=xinp+xinp.*randn(size(xinp)).*data.cnoise;
            else
%                 xinp  =(([predc(t-1,:) data.age_h(t-1)]-data.xmin)./(data.xmax-data.xmin))*2-1;%João Ramos
                xinp  =(([predc(t-1,:)]-data.xmin)./(data.xmax-data.xmin))*2-1;%João Ramos
            end
            hnet=hnetff(hnet,xinp);
            
            if hnet.ymode~=1  %João Ramos
               mx(t-1)= predc(t-1,1)*data.deltat*data.batch(ib).vol(t-1)/1000;
            end
            
            mr = mr + mx(t-1)*(stoichiometry*hnet.y);
            dmrdy = dmrdy + mx(t-1)*stoichiometry;
            if hnet.ymode~=1  %João Ramos
                predc(t,:)=(data.batch(ib).accum(t,:)+mr')/(data.batch(ib).vol(t)*1e-3);
            end
            
            if rand(1,1) <= hnet.minibatch
                
                count = count+numel(mr);
                ytarg = ybatch(t,:)';
                if hnet.ymode==1  %João Ramos
                    res  = (ytarg-mr)./data.ystd(:,ib);
                    res2 = (ytarg-mr)./data.ystd2(:,ib);
                else
                    res  = (ytarg- predc(t,:)')./data.ystd(:,ib); %João Ramos
                    res2 =(ytarg- predc(t,:)')./data.ystd2(:,ib);
                end
                dlossdy = (-2*res./data.ystd(:,ib))'*dmrdy;
                mse = mse + res'*res ;
                mse2= mse2 + res2'*res2 ;
                
                % penalty for negative metabolic fluxes
                if strcmp(hnet.layers{hnet.nl}.type,'Metabolic') 
                    vall = hnet.layers{hnet.nl}.v;
                    resv = min(vall(ir),0);
%                     vir_=max(vall(ir),0);
%                     vall(ir)=vir_;
%                     accumin=siir*vall;
%                     resv = log(cosh(vir));
                    lambda=1e6;
                    mse = mse + lambda*(resv')*resv;
                    hnet.layers{hnet.nl}.dlossdx = ...
                        hnet.layers{hnet.nl}.dlossdx +nsiirt*(2*lambda*resv);
                end
                hnet=hnetbp(hnet,dlossdy');
            end
            
        end
    end
end
if mse==0 % no batch selected due to minibatch
    mse=inf;
end
end


function hnet_verify_grads(hnet,data) %#####################################
%minibatch = hnet.minibatch;
%dropout = hnet.dropout;
hnet.minibatch=1;
hnet.dropout=0;
fobj=@(w)lossfunvec(w,hnet,data);
winit=hnetgetw(hnet);

options = optimoptions(@fminunc,'Algorithm','trust-region',...
    'CheckGradients',true,'SpecifyObjectiveGradient',true,...
    'Display','iter','MaxIter',1);  
[wfinal,fval,res,exitflag] = fminunc(fobj,winit,options);

end

function [loss,grads]=lossfunvec(w,hnet,data)%#############################

hnet=hnetsetw(hnet,w);
[hnet,loss,count]=lossfun(hnet,data);
[~,sum_dlossdw]=hnetgetgrad(hnet);
grads = sum_dlossdw';

end
