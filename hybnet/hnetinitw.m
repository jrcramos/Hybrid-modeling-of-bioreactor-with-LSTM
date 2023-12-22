function hnet = hnetinitw(hnet)%###########################################
for i=1:hnet.nl
    hnet.layers{i}=feval(hnet.layers{i}.initwfun, hnet.layers{i});
end
end