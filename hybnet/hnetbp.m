function hnet=hnetbp(hnet,dlossdy)%################################
% hnet.dlossdy =  dlossdy;
% for i=1:hnet.nl
%     hnet.layers{i}.dlossdy=zeros(hnet.layers{i}.ny,1);
% end
hnet.layers{hnet.nl}.dlossdy = dlossdy;
for i = hnet.nl:-1:1
    hnet.layers{i}=feval(hnet.layers{i}.bpfun,hnet.layers{i},hnet.layers{i}.dlossdy);
    inpl=hnet.layers{i}.inpl;
    if inpl==0
        %hnet.dlossdx = hnet.dlossdx + hnet.layers{i}.dlossdx;%???
    else
        hnet.layers{inpl}.dlossdy = hnet.layers{inpl}.dlossdy + ...
            hnet.layers{i}.dlossdx;
    end
end
end
