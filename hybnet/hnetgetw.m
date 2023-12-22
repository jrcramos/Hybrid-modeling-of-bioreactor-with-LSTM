function [w]=hnetgetw(hnet)%###############################################
w=zeros(hnet.nw,1);
count=1;
for i=1:hnet.nl
    w(count:count+hnet.layers{i}.nw-1) = ...
        feval(hnet.layers{i}.getwfun,hnet.layers{i});
    count = count + hnet.layers{i}.nw;
end
end