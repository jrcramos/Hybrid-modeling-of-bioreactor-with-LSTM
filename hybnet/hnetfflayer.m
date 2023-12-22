function layer=hnetfflayer(ninputs,numnodes,node,inpl)%########################
if nargin<4
    inpl=-1;
end

layer.inpl = inpl;
layer.nx = ninputs;
layer.ny = numnodes;
layer.nw = layer.ny*layer.nx+layer.ny;
%initialize gradients
layer.dlossdw=zeros(layer.ny,layer.nx);
layer.dlossdb=zeros(layer.ny,1);
layer.sum_dlossdw=zeros(layer.ny,layer.nx);
layer.sum_dlossdb=zeros(layer.ny,1);
%initialize adam
layer.m_w=zeros(layer.ny,layer.nx);
layer.m_b=zeros(layer.ny,1);
layer.v_w=zeros(layer.ny,layer.nx);
layer.v_b=zeros(layer.ny,1);
%???????????????????
%layer = structfun(@gpuArray, layer, 'UniformOutput', false);

layer.type = 'feedforward';
layer.node = node;
layer.initwfun = @ffinitw;
layer.setwfun = @ffsetw;
layer.getwfun = @ffgetw;
layer.getgradfun = @ffgetgrad;
layer.fffun = @ffff;
layer.bpfun = @ffbp;
layer.bpresetfun = @ffresetbp;
layer.updwfun=@ffupdweights;
layer.initstatefun = @ffinitstate;

[layer.nodefffun,layer.nodebpfun]=hnetnode(node);
%initialize weights
dw = 2; %%%1 Normal  %%%2 Xavier
layer = ffinitw(layer,dw);

end

function layer=ffinitw(layer,dw)%########################################
if nargin<2 
    dw=1;
end
if dw==1
    dw=0.25;
    %%1 normal
    layer.w = randn(layer.ny,layer.nx)*dw;
    layer.b = randn(layer.ny,1)*dw;
else
    %%%2 Xavier
    n=layer.nw;%layer.ny*layer.nx+layer.ny;
    lower = -(1.0 / sqrt(n));
    upper =  (1.0 / sqrt(n));
    layer.w = lower + rand(layer.ny,layer.nx) * (upper - lower);
    layer.b = lower + rand(layer.ny,1) * (upper - lower);
end

end

function layer=ffinitstate(layer)%########################################
end

function layer=ffsetw(layer,w)%############################################
count=layer.ny*layer.nx;
layer.w = reshape(w(1:count),layer.ny,layer.nx);
layer.b = reshape(w(count+1:count+layer.ny),layer.ny,1);
end

function w=ffgetw(layer)%##################################################
count = layer.ny*layer.nx;
w = zeros(layer.nw,1);
w(1:count) = reshape(layer.w,count,1);
w(count+1:count+layer.ny) = reshape(layer.b,layer.ny,1);
end

function [dlossdw,sum_dlossdw]=ffgetgrad(layer)%###############################
dlossdw=zeros(layer.nw,1);
count=layer.ny*layer.nx;
dlossdw(1:count) = reshape(layer.dlossdw,count,1);
dlossdw(count+1:count+layer.ny) = reshape(layer.dlossdb,layer.ny,1);

sum_dlossdw=zeros(layer.nw,1);
count=layer.ny*layer.nx;
sum_dlossdw(1:count) = reshape(layer.sum_dlossdw,count,1);
sum_dlossdw(count+1:count+layer.ny) = reshape(layer.sum_dlossdb,layer.ny,1);
end

function layer=ffff(layer,x)%##############################################
layer.x=x;
layer.xh=layer.w*layer.x+layer.b;
layer.y=feval(layer.nodefffun,layer.xh);
end

function layer=ffbp(layer,dlossdy)%########################################
dydxh=feval(layer.nodebpfun,layer.xh,layer.y);
dlossdxh=dlossdy.*dydxh;
layer.dlossdw=dlossdxh*layer.x';
layer.dlossdb=dlossdxh;
layer.dlossdx=layer.w'*dlossdxh;
layer.sum_dlossdw = layer.sum_dlossdw + layer.dlossdw;
layer.sum_dlossdb = layer.sum_dlossdb + layer.dlossdb;
end

function layer=ffresetbp(layer)%########################################
layer.sum_dlossdw=zeros(layer.ny,layer.nx);
layer.sum_dlossdb=zeros(layer.ny,1);
end

function layer=ffupdweights(layer,iter,alfa,beta1,beta2,eta,dropout) %#####
%adam update rules

if iter <=1 %initialize adam moments
    layer.m_w=zeros(layer.ny,layer.nx);
    layer.m_b=zeros(layer.ny,1);
    layer.v_w=zeros(layer.ny,layer.nx);
    layer.v_b=zeros(layer.ny,1);
end

ind = rand(layer.ny,1)>=dropout;  %nodes to be updated


layer.m_w(ind,:) = beta1*layer.m_w(ind,:) + (1-beta1)*layer.sum_dlossdw(ind,:);
layer.m_b(ind)  =  beta1*layer.m_b(ind)   + (1-beta1)*layer.sum_dlossdb(ind);

layer.v_w(ind,:) = beta2*layer.v_w(ind,:) + (1-beta2)*(layer.sum_dlossdw(ind,:).*layer.sum_dlossdw(ind,:));
layer.v_b(ind)   = beta2*layer.v_b(ind)   + (1-beta2)*(layer.sum_dlossdb(ind).*layer.sum_dlossdb(ind));

mest_w = layer.m_w(ind,:)/(1-beta1^iter);
mest_b = layer.m_b(ind)/(1-beta1^iter);

vest_w = layer.v_w(ind,:)/(1-beta2^iter);
vest_b = layer.v_b(ind)/(1-beta2^iter);

layer.w(ind,:) = layer.w(ind,:) - alfa*(mest_w./(sqrt(vest_w)+eta));
layer.b(ind) = layer.b(ind) - alfa*(mest_b./(sqrt(vest_b)+eta));

end