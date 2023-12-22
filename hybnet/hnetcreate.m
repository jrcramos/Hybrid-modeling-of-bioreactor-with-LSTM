function hnet=hnetcreate(layers,tag)%######################################
if nargin<2
    tag='HNET';
end
hnet.tag = tag;
hnet.layers = layers;
hnet.nl = numel(layers);
hnet.nx = layers{1}.nx;
hnet.ny = layers{hnet.nl}.ny;
hnet.nw = 0;
for i=1:hnet.nl
    hnet.nw = hnet.nw + hnet.layers{i}.nw;
end
hnet.setwfun=@hnetsetw;
hnet.getgradfun=@hnetgetgrad;
hnet.fffun=@hnetff;
hnet.bpfun=@hnetbp;
%adam updae rules
hnet.iter=0;
hnet.alfa = 0.001; %AUTOMATE, if undefined, make it n
hnet.beta1=0.9;
hnet.beta2=0.999;
hnet.eta=1e-8;
hnet.minibatch = 0.8;
hnet.dropout = 0.1;
%verify connectivity
for i=1:hnet.nl
    if hnet.layers{i}.inpl==-1
        hnet.layers{i}.inpl=i-1;
    else
        ind = hnet.layers{i}.inpl>=i;
        assert(sum(ind)==0,'layer connectivy error')
    end
end
end








