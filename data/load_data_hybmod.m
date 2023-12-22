function [x,y,accum,formed,vol,age_h,xlabels,ylabels,reacted_masses]=load_data_hybmod(ymode,varargin)
load data.mat
if nargin>1
    load([varargin{1} '.mat']);
end

age_h=data(1).time;

names_met=lower(data(1).names);


xlabels=string(names_met);
ylabels=string(names_met(1:numel(names_met)));

x={};y={};reacted_masses=[];accum={};vol={};formed=[];
for i=1:size(data,2)
    x=[x {[data(i).conc]}];
    formed=[formed;{data(i).m_r}];
    reacted_masses=[reacted_masses; data(i).m_r'];
    if ymode==1
        y=[y {data(i).m_r}];
    else
        y=[y {data(i).conc}];
    end
    accum=[accum {data(i).accum'}]; 
    vol=[vol {data(i).vol}];
end
end
