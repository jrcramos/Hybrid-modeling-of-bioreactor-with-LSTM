function hnet=hnet_update_weights(hnet)%################################
hnet.iter = hnet.iter+1;
hnet.alfa=hnet.alfa0-(hnet.alfa0-hnet.alfaf)*hnet.iter/hnet.nter;%hnet.alfa0./sqrt(hnet.iter/50);% Joao Ramos
for i = 1:hnet.nl
    if hnet.layers{i}.nw > 0
        hnet.layers{i}=feval(hnet.layers{i}.updwfun,hnet.layers{i},...
            hnet.iter,hnet.alfa,hnet.beta1,hnet.beta2,hnet.eta,hnet.dropout);
    end
end
end
