
function [batch]=ReadBatchData()
batch=struct(); batch.nbatch=9;

%---------------------------------------------------------------------------
% first preprocessing, file names and sheet names
%--------------------------------------------------------------------------
filenames='data.xlsx';
sheets_data1='synthetic data';% sheets the get data from
% sheets_data2=sheetnames(filenames{2});sheets_data2=cellstr(sheets_data2(2:end));% sheets the get data from


%---------------------------------------------------------------------------
% Get a single time period for all batchs and all data
%--------------------------------------------------------------------------
fprintf('___________________________________________________________________________\n\n')
fprintf('Building a unique time vector \n', 0)
[num,txt,raw]=xlsread(filenames,sheets_data1);

time=unique([0; cell2mat(raw(3:end,1))]);
batch.age=time;
br="";
for i=1:9
    br(i,1)="Br" + string(num2str(i));
end
batch.batchid=br;
line_col=[find(cell2mat(raw(3:end,1))==0)+2 find(cell2mat(raw(3:end,1))==240)+2];

%---------------------------------------------------------------------------
% Get VCD,mAb, Glc, aa, and feed volumes
%--------------------------------------------------------------------------
fprintf('___________________________________________________________________________\n')
fprintf('Getting conc and feeds from %s %2.2f%% concluded.\n', filenames, 0)
for i=2:31 %
   
    temp=NaN(numel(time),9); % creat a matrix sized time x nbatch
    for ib=1:9  % find data each batch
        temp(:,ib)=cell2mat(raw(line_col(ib,1):line_col(ib,2),i));
    end
    if i==3 % save product
        batch.product(1).val=temp;
        batch.product(1).unit={'mg'};%raw(2,1);
    elseif i==2 % save VCD
        batch.vcd(1).val=temp;
        batch.vcd(1).unit={'Mcell/mL'};
    elseif i>=3 && i<=26 % save conc
        batch.met(i-3).conc=temp;
        batch.met(i-3).name=raw(1,i);
        batch.met(i-3).unit=raw(2,i);
    elseif i==30
        batch.samples.samplevol=temp; % non cumulative sample
    elseif i==31  % save feed data
        batch.feed.feedvol=cumsum(temp);
        batch.feed.feedid='feed1';
    end
    clc
    fprintf('Getting conc and feeds from %s %2.2f%% concluded.\n', filenames, i/31*100)
end


%---------------------------------------------------------------------------
% Estimate bioreactor volume from feed and sample+bleed volumes
%--------------------------------------------------------------------------
fprintf('Estimating bioreactor volume based on sampling+bleed and feeding.\n')
temp=NaN(numel(time),9); % creat a matrix sized time x nbatch
for ib=1:9  % find data each batch
    temp(:,ib)=cell2mat(raw(line_col(ib,1):line_col(ib,2),27));
end
batch.samples.reactorvolume=temp;



%---------------------------------------------------------------------------
% Get feed recepies and match to feed volumes previously imported
%--------------------------------------------------------------------------
fprintf('Getting feed recepies from %s.\n', filenames)
[num,txt,raw]=xlsread(filenames,'synthetic DoE');

for i=1
    metnames=deblank(raw(12,8:30));
    unit=deblank(raw(13,8:30));
    vals=cell2mat(raw(14,8:30));
    batch.feedcomp(i).id='feed1';
    batch.feed(i).recepyname='feed1';
    batch.feedcomp(i).conc=vals;
    batch.feedcomp(i).mets=metnames;
    batch.feedcomp(i).unit=unit;
end


end


