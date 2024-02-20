function [batch]=CalcReactedAmount(batch)
nbatch=batch.nbatch;

%--------------------------------------------------------------------------
% Built time vector per batch where conc, feed, sampling occurrs
%--------------------------------------------------------------------------
for i=1:nbatch
    datasbatch=batch.vcd.val(:,i);datasbatch(:,2)=batch.product.val(:,i);datasbatch(:,3)=batch.samples.samplevol(:,i);datasbatch(:,4)=batch.met(7).conc(:,i);
    datasbatch(:,5)=batch.feed(1).feedvol(:,i);
    for ib=1:size(batch.feed,2)
        datasbatch(:,end+1)=batch.feed(ib).feedvol(:,i);
    end
    datasbatch(datasbatch<=1e-6)=NaN; datasbatch=nansum(datasbatch,2);
    idst=find(datasbatch>0);
    
    batch.perbatch(i).age=batch.age(idst);
    batch.perbatch(i).ids=idst;
    
end



%--------------------------------------------------------------------------
% Compute ICV - Integral of Viable Cell density
%--------------------------------------------------------------------------
batch.age_all=linspace(0,batch.age(end),150);

for i=1:nbatch
   % get measured VCD, vol 
   ids= batch.perbatch(i).ids;
   vcdvol=batch.vcd.val(ids,i).*batch.samples.reactorvolume(ids,i);
   age = batch.perbatch(i).age;
   ind = ~isnan(vcdvol);   
   vcdvol=vcdvol(ind);  
   age=age(ind);
   
   % calculate IVC    
   batch.ivc(i).val = zeros(length(batch.perbatch(i).age),1);
   logxv=interp1(age,log(vcdvol),batch.perbatch(i).age,'linear','extrap');  
   xv = exp(logxv); %fill the gaps
   for j=2:length(batch.perbatch(i).age)
       dt=batch.perbatch(i).age(j)-batch.perbatch(i).age(j-1);
       if dt==0
           fprintf('Fatal error!!!! time points with the same age\n')
           return
       end
       mu=log(xv(j)/xv(j-1))/dt;
       if abs(mu)>0
            xvmed=(xv(j)-xv(j-1))/mu/dt;
       else
            xvmed=xv(j);
       end
       batch.ivc(i).val(j)=batch.ivc(i).val(j-1) + xvmed*dt;
   end

   % calculate IVC_all 
   batch.ivc_all(i).val = zeros(length(batch.age_all),1);
   logxv=interp1(age,log(vcdvol),batch.age_all,'linear','extrap');  
   xv = exp(logxv);
   for j=2:length(batch.age_all)
       dt=batch.age_all(j)-batch.age_all(j-1);
       if dt==0
           fprintf('Fatal error!!!! time points with the same age\n')
           return
       end
       mu=log(xv(j)/xv(j-1))/dt;
       if abs(mu)>0
            xvmed=(xv(j)-xv(j-1))/mu/dt;
       else
            xvmed=xv(j);
       end
       batch.ivc_all(i).val(j)=batch.ivc_all(i).val(j-1) + xvmed*dt;
   end
   
end



%--------------------------------------------------------------------------
% Compute amount of mAb and metabolite IN to the reactor
%--------------------------------------------------------------------------
% replace NaN samples with 0 and correct feed
batch.samples.samplevol(isnan(batch.samples.samplevol))=0;
batch.samples.samplevol=cumsum(batch.samples.samplevol);
for ib=1:size(batch.feed,2)
   batch.feed(ib).feedvol(1,:)=0;
   batch.feed(ib).feedvol=fillmissing(batch.feed(ib).feedvol,'previous');
end
%%%

fprintf('___________________________________________________________________________\n')
fprintf('Compute amount of mAb and metabolite IN to the reactor %2.2f%% concluded.\n', 0)
for i=1:nbatch
    ids= batch.perbatch(i).ids;
    n =length(batch.perbatch(i).age);n2=numel(batch.age);
    batch.vcd_in.val(:,i)=NaN(n2,1);
    batch.product_in.val(:,i)=NaN(n2,1);
    batch.vcd_in.val(1,i) =batch.vcd.val(1,i).*batch.samples.reactorvolume(1,i); %inicial amount
    batch.product_in.val(1,i) = batch.product.val(i,1)*batch.samples.reactorvolume(1,i)*1-3;
    batch.product_in.val(1,i)=0; %???????????????   product has nan in the begining

    for j=2:n   %only sample event influencies VCD and mAb, i.e. not in the feeds ever
        vcd=batch.vcd.val(ids(j),i);
        if isnan(vcd)
            vcd=0; %problem missing vaue
        end
        sample_v=batch.samples.samplevol(ids(j),i)-max(batch.samples.samplevol(1:ids(j-1),i));
              
        batch.vcd_in.val(ids(j),i) = batch.vcd_in.val(ids(j-1),i) - ...
            sample_v*vcd;
        product=batch.product.val(ids(j),i);
        if isnan(product)
            product=0;  %problem missing value
        end
        batch.product_in.val(ids(j),i) = batch.product_in.val(ids(j-1),i) - ...
            sample_v*1e-3*product;
    end
    
    
     for k=1:size(batch.met,2)  %metabolites
         
         batch.met_in(k).val(:,i)=NaN(n2,1);
          
         if ~isnan(batch.met(k).conc(1,i))
             batch.met_in(k).val(1,i) = batch.met(k).conc(1,i)*batch.samples.reactorvolume(1,i)*1e-3; %initial amount
         else
             fprintf('Error NaN at first position \n')
             pause
         end
         
          for j=2:n
              feed_data=0;
              for ifeed=1:size(batch.feedcomp,2)
                  if find(ismember(batch.met(k).name,batch.feedcomp(ifeed).mets))
                      idm=find(ismember(batch.feedcomp(ifeed).mets,batch.met(k).name));
                      fconc=batch.feedcomp(ifeed).conc(idm);
                      idf=find(strcmp(batch.feedcomp(ifeed).id,[batch.feed.recepyname]));
                      if ~isempty(idf) && batch.feed(idf).feedvol(ids(j),i)>0
                          fedvol=batch.feed(idf).feedvol(ids(j),i)-max(batch.feed(idf).feedvol(1:ids(j-1),i));% because feed volume is cumulative
                          feed_data=feed_data+fedvol*1e-3*fconc;
                      end
                  end
              end
              sample_v=batch.samples.samplevol(ids(j),i)-max(batch.samples.samplevol(1:ids(j-1),i));
              idn= find(~isnan(batch.met(k).conc(:,i)));
              [~,idc]=min(batch.age(idn)-batch.age(ids(j)));
              sample_data=sample_v*1e-3*batch.met(k).conc(idn(idc),i); % just remove mass based in closest measurement
              
                      
              batch.met_in(k).val(ids(j),i) =  batch.met_in(k).val(ids(j-1),i) ...
                        -sample_data  + feed_data;
          end
          batch.met_in(k).name=batch.met(k).name;  
     end
     disp([repmat(char(8), 1, 75)])
     fprintf('Compute amount of mAb and metabolite IN to the reactor %2.2f%% concluded.\n', i/nbatch*100)
end





%--------------------------------------------------------------------------
% Computate amount of species that reacted
%--------------------------------------------------------------------------
for i=1:nbatch
   batch.m_r_vcd(:,i)     = batch.vcd.val(:,i).*batch.samples.reactorvolume(:,i) - batch.vcd_in.val(:,i);
   batch.m_r_product(:,i) = batch.product.val(:,i).*batch.samples.reactorvolume(:,i)*1e-3 - batch.product_in.val(:,i);   
   for k=1:size(batch.met,2)
%        if isnan(batch.met(k).conc(1,i))
%            batch.met(k).conc(1,i)=1.0*batch.met(k).conc(2,nbatch); % initial conc was measured in batch nbatch
%        end
       batch.met_m_r(k).val(:,i) = batch.met(k).conc(:,i).*batch.samples.reactorvolume(:,i)*1e-3 - batch.met_in(k).val(:,i);
       batch.met_m_r(k).name=batch.met(k).name; 
   end
end





end



