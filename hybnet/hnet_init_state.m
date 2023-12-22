function hnet=hnet_init_state(hnet)%################################
for i = 1:hnet.nl  
  hnet.layers{i}=feval(hnet.layers{i}.initstatefun,hnet.layers{i});
end
end
