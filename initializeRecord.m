function initializeRecord()
%INITIALIZERECORD Summary of this function goes here
%   Detailed explanation goes here
global groupbuffer;
global idlist;
global msglist;
global fsynctype;
global groupID;
global groupTOA;
global bWrt;
global LIGHTSPEED;
LIGHTSPEED = physconst('LightSpeed')*1e-3;
global jd2000;
date2000=datetime('1-Jan-2000 12:00:00');
jd2000=juliandate(date2000);

groupbuffer=cell(10,25,100);
idlist=cell(1,100);
msglist=cell(1,100);
fsynctype = zeros(1,100);
groupID=zeros(1,100);
groupTOA=repmat(datetime,1,100);
bWrt=zeros(1,100);
end

