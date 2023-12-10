function readrinex(filename)
%READRINEX Summary of this function goes here
%   Detailed explanation goes here
global ephdata;
ephdata=cell(1,36);
info = rinexread(filename);
data=info.Galileo;
for i=1:36
idx=data.SatelliteID==i;
ephdata{i}=data(idx,:);
end
end

