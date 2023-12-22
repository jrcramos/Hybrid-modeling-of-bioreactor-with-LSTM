function [dlossdw,sum_dlossdw]=hnetgetgrad(hnet)%##############################
dlossdw    =zeros(hnet.nw,1);
sum_dlossdw=zeros(hnet.nw,1);
count=1;
for i=1:hnet.nl
    [dlossdw_i,sum_dlossdw_i] = feval(hnet.layers{i}.getgradfun,hnet.layers{i});
    dlossdw(count:count+hnet.layers{i}.nw-1)    = dlossdw_i;
    sum_dlossdw(count:count+hnet.layers{i}.nw-1)= sum_dlossdw_i;    
    count = count + hnet.layers{i}.nw;
end
end