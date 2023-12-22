function hnet=hnet_reset_bp(hnet)%################################
for i = 1:hnet.nl  
  hnet.layers{i}=feval(hnet.layers{i}.bpresetfun,hnet.layers{i});
end
end
