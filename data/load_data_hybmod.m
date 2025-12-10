function [x,y,accum,formed,vol,age_h,xlabels,ylabels,reacted_masses,s_extMatch,s_int,isrev]=load_data_hybmod(ymode,varargin)
load data.mat
load model.mat
if nargin>1
    load([varargin{1} '.mat']);
end


Aint=model.Sint;
Aext=model.Sext;
metstr_ext = model.mets(~model.isintr);
mets_data=data(1).names;
ids=[];
for i=1:numel(mets_data)
    ids=[ids strmatch(mets_data{i},metstr_ext)];
end
s_extMatch=Aext(ids,:);
s_int=Aint;
isrev=model.isrev;


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
