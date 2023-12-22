function hnet=hnetff(hnet,inpdata)%########################################
inpdata = reshape(inpdata,numel(inpdata),1);
hnet.x = inpdata;
for i=1:hnet.nl
    inpl=hnet.layers{i}.inpl;
    if inpl==0
        inp = inpdata;
    else
        inp = hnet.layers{inpl}.y;
    end
    hnet.layers{i}=feval(hnet.layers{i}.fffun,hnet.layers{i},inp);
end
hnet.y=hnet.layers{hnet.nl}.y;

% revise in the future
hnet.dlossdy =  zeros(hnet.ny,1);
for i=1:hnet.nl
    hnet.layers{i}.dlossdy=zeros(hnet.layers{i}.ny,1);
end

end
