function layer=hnetLSTMlayer(ninputs,numnodes,inpl)%#######################
if nargin<3
    inpl=-1;
end
layer.type = 'LSTM';
layer.nx = ninputs;
layer.ny = numnodes;
layer.inpl = inpl; 
layer.y=zeros(numnodes,1);
layer.cold = zeros(numnodes,1);
layer.dcdf= zeros(numnodes,1);
layer.dcdi= zeros(numnodes,1);
layer.dcdcest= zeros(numnodes,1);
layer.c=zeros(layer.ny,1);

layer.fg = hnetfflayer(ninputs+numnodes,numnodes,'sigm');
layer.ig = hnetfflayer(ninputs+numnodes,numnodes,'sigm');
layer.og = hnetfflayer(ninputs+numnodes,numnodes,'sigm');
layer.ci = hnetfflayer(ninputs,numnodes,'tanh');

layer.nw = layer.fg.nw+layer.ig.nw+layer.og.nw+layer.ci.nw;
layer.initwfun = @lstminitw;
layer.setwfun = @lstmsetw;
layer.getwfun = @lstmgetw;
layer.getgradfun = @lstmgetgrad;
layer.fffun = @lstmff;
layer.bpfun = @lstmbp;
layer.bpresetfun = @lstmresetbp;
layer.updwfun=@lstmupdweights;
layer.initstatefun=@lstminitstate;
end

function layer=lstminitw(layer,dw)%########################################
if nargin<2
    dw=2; %%%1 Normal  %%%2 Xavier
end

layer.fg = feval(layer.fg.initwfun,layer.fg,dw);
layer.ig = feval(layer.ig.initwfun,layer.ig,dw);
layer.og = feval(layer.og.initwfun,layer.og,dw);
layer.ci = feval(layer.ci.initwfun,layer.ci,dw);

end

function layer=lstminitstate(layer)%########################################
layer.c=zeros(layer.ny,1);
layer.cold=zeros(layer.ny,1);
layer.y=zeros(layer.ny,1);
layer.dcdf= zeros(layer.ny,1);
layer.dcdi= zeros(layer.ny,1);
layer.dcdcest= zeros(layer.ny,1);
end

function layer=lstmsetw(layer,w)%############################################
count = 1;
layer.fg = feval(layer.fg.setwfun,layer.fg,w(count:count+layer.fg.nw-1,1));
count = count + layer.fg.nw;
layer.ig = feval(layer.ig.setwfun,layer.ig,w(count:count+layer.ig.nw-1,1));
count = count + layer.ig.nw;
layer.og = feval(layer.og.setwfun,layer.og,w(count:count+layer.og.nw-1,1));
count = count + layer.og.nw;
layer.ci = feval(layer.ci.setwfun,layer.ci,w(count:count+layer.ci.nw-1,1));
end

function w=lstmgetw(layer)%################################################
w = zeros(layer.nw,1);
count = 1;
w(count:count+layer.fg.nw-1,1) = feval(layer.fg.getwfun,layer.fg);
count = count + layer.fg.nw;
w(count:count+layer.ig.nw-1,1) = feval(layer.ig.getwfun,layer.ig);
count = count + layer.ig.nw;
w(count:count+layer.og.nw-1,1) = feval(layer.og.getwfun,layer.og);
count = count + layer.og.nw;
w(count:count+layer.ci.nw-1,1) = feval(layer.ci.getwfun,layer.ci);
end

function [dlossdw,sum_dlossdw]=lstmgetgrad(layer)%#########################
dlossdw=zeros(layer.nw,1);
sum_dlossdw=zeros(layer.nw,1);
count = 1;
[dlossdw(count:count+layer.fg.nw-1,1),...
    sum_dlossdw(count:count+layer.fg.nw-1,1)] = ...
    feval(layer.fg.getgradfun,layer.fg);
count = count + layer.fg.nw;
[dlossdw(count:count+layer.ig.nw-1,1),...
    sum_dlossdw(count:count+layer.ig.nw-1,1)] = ...
    feval(layer.ig.getgradfun,layer.ig);
count = count + layer.ig.nw;
[dlossdw(count:count+layer.og.nw-1,1),...
    sum_dlossdw(count:count+layer.og.nw-1,1)] = ...
    feval(layer.og.getgradfun,layer.og);
count = count + layer.og.nw;
[dlossdw(count:count+layer.ci.nw-1,1),...
    sum_dlossdw(count:count+layer.ci.nw-1,1)] = ...
    feval(layer.ci.getgradfun,layer.ci);
end

function layer=lstmff(layer,x)%############################################
layer.cold=layer.c;  %??????

layer.x=x;

inp = [x;layer.cold];

layer.fg = feval(layer.fg.fffun,layer.fg,inp);
layer.f = layer.fg.y;
layer.ig = feval(layer.ig.fffun,layer.ig,inp);
layer.i=layer.ig.y;
layer.og = feval(layer.og.fffun,layer.og,inp);
layer.o = layer.og.y; 
layer.ci = feval(layer.ci.fffun,layer.ci,x);
layer.cest = layer.ci.y;

layer.c= layer.f.*layer.cold + layer.i.*layer.cest;

layer.y = layer.o.*layer.c;

layer.dcdf=      layer.f.*layer.dcdf    +   layer.cold;
layer.dcdi=      layer.f.*layer.dcdi    +   layer.cest;
layer.dcdcest=   layer.f.*layer.dcdcest +   layer.i;

% layer.dcdf=      layer.cold;
% layer.dcdi=      layer.cest;
% layer.dcdcest=   layer.i;

end

function layer=lstmbp(layer,dlossdy)%########################################
%y = o.*c; ------
dlossdo = dlossdy.*layer.c;
dlossdc = dlossdy.*layer.o;
%c= f.*cold + i.*cest; ----------
dlossdf = dlossdc.*layer.dcdf;
dlossdi = dlossdc.*layer.dcdi;
dlossdcest = dlossdc.*layer.dcdcest;

%layer.ci = feval(layer.ci.fffun,layer.ci,inp);
layer.ci = feval(layer.ci.bpfun,layer.ci,dlossdcest);
%layer.og = feval(layer.og.fffun,layer.og,inp);
layer.og = feval(layer.og.bpfun,layer.og,dlossdo);
%layer.ig = feval(layer.ig.fffun,layer.ig,inp);
layer.ig = feval(layer.ig.bpfun,layer.ig,dlossdi);
%layer.fg = feval(layer.fg.fffun,layer.fg,inp);
layer.fg = feval(layer.fg.bpfun,layer.fg,dlossdf);
layer.dlossdx=layer.ci.dlossdx(1:layer.nx,1)+...
    layer.og.dlossdx(1:layer.nx,1)+...
    layer.ig.dlossdx(1:layer.nx,1)+...
    layer.fg.dlossdx(1:layer.nx,1);
layer.dlossdcold=layer.og.dlossdx(1+layer.nx:end,1)+...
    layer.ig.dlossdx(1+layer.nx:end,1)+...
    layer.fg.dlossdx(1+layer.nx:end,1);
end

function layer=lstmresetbp(layer)%########################################
layer.ci=feval(layer.ci.bpresetfun,layer.ci);
layer.og=feval(layer.og.bpresetfun,layer.og);
layer.ig=feval(layer.ci.bpresetfun,layer.ig);
layer.fg=feval(layer.fg.bpresetfun,layer.fg);
end

function layer=lstmupdweights(layer,iter,alfa,beta1,beta2,eta,dropout)%#####
%adam update rules
layer.fg=feval(layer.fg.updwfun,layer.fg,iter,alfa,beta1,...
    beta2,eta,dropout);
layer.ig=feval(layer.ig.updwfun,layer.ig,iter,alfa,beta1,...
    beta2,eta,dropout);
layer.og=feval(layer.og.updwfun,layer.og,iter,alfa,beta1,...
    beta2,eta,dropout);
layer.ci=feval(layer.ci.updwfun,layer.ci,iter,alfa,beta1,...
    beta2,eta,dropout);
end