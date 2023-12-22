function [fffun,bpfun]=hnetnode(node)
switch node
    case 'tanh'
        fffun=@tanhff;
        bpfun=@tanhbp;
    case 'sigm'
        fffun=@sigmff;
        bpfun=@sigmbp;
    case 'relu'
        fffun=@reluff;
        bpfun=@relubp;
    case 'linear'
        fffun=@linearff;
        bpfun=@linearbp;
end
end


%tanh activation node-------------------------------------------------
function y=tanhff(x)
y=tanh(x);
end
function yl=tanhbp(x,y)
yl=1-y.*y;
end


%sigmoid activation function-------------------------------------------------
function y=sigmff(x)
y=1./(1+exp(-x));
end
function yl=sigmbp(x,y)
yl=y.*(1-y);
end

%relu activation function-------------------------------------------------
function y=reluff(x)
y=max(0.01*x,0.003*x);
end

function yl=relubp(x,y)
yl=zeros(size(y));
yl(y>=0)=0.01;
yl(y<0)=0.003;
end

%linear activation function-----------------------------------------------
function y=linearff(x)
y=x;
end
function yl=linearbp(x,y)
yl=ones(size(y));
end
