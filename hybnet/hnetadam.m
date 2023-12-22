function [w,fval]=hnetadam(lossfunc, w0, options, ofun)
% default adam parameters; in future goes to structure options
if nargin<3
    options=[];
end
if isfield(options,"adalfa")
    alfa=  options.adalfa;
else
    alfa = 0.001; %AUTOMATE, if undefined, make it n
end
beta1=0.9;
beta2=0.999;
eta=1e-8;
if isfield(options,"admini")
    minibatchsize=options.admini;
else
    minibatchsize = 0.8; %AUTOMATE, if undefined, make it n
end
if isfield(options,"addrop")
    dropout=options.addrop;
else
    dropout = 0.2; %AUTOMATE, if undefined, make it n
end
%initialization of variables
w=reshape(w0,numel(w0),1); %weights
m=zeros(numel(w),1);  %first moment vector
v=zeros(numel(w),1); %Second moment vector


for i=1:options.niter  %iteration cylce
  
   witer=w;
   
  %call lossfun get residuals and gradients residuals
  [fval,g]=feval(lossfunc,witer,minibatchsize,dropout);  

   m = beta1*m + (1-beta1)*g;
   v = beta2*v + (1-beta2)*(g.*g);  
   mest = m/(1-beta1^i);
   vest = v/(1-beta2^i);
   w = w - alfa*(mest./(sqrt(vest)+eta));
  
end
[fval]=feval(lossfunc,witer,minibatchsize,dropout);  


end

