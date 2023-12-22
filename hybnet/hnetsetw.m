function hnet=hnetsetw(hnet,w)%###########################################
count=1;
for i=1:hnet.nl
    hnet.layers{i}=feval(hnet.layers{i}.setwfun, hnet.layers{i},...
        w(count:count+hnet.layers{i}.nw-1));
    count = count + hnet.layers{i}.nw;
end
end