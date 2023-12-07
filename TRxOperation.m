function [TOA,FOA] = TRxOperation(satIDs,ToT,FoT,TxSite,RxSite)

validChns=find(satIDs);
satIDs=satIDs(validChns);
noOfSats=length(validChns);
t0=split2fields(repmat(ToT,1,noOfSats));
f_list=[1544.1e6,1544.9e6,1544.21e6];
%uplink
t1=addSeconds(t0,0.16);
[posS,velS,dt1]=actualtof(t1,satIDs,TxSite,'uplink');
t2=addSeconds(t0,0.08);
[posS1,velS1,~]=actualtof(t2,satIDs,TxSite,'uplink');%for doppler calculation
fd1=getDoppler(posS1,velS1,TxSite,FoT);
%downlink
dt2=tof(posS,RxSite);
cflag=floor(satIDs/100)-3;
freq_trns=f_list(cflag) - 406.05e6;
fc1=FoT+fd1+freq_trns;
fd2 = getDoppler(posS1,velS1,RxSite,fc1);
%total with noise added
TOA=ToT+(dt1+dt2+20e-6*randn(1,noOfSats))/86400;
FOA=fc1+fd2;
end

function[posS,velS,dt]= actualtof(t,sids,place,journey)
dt=0;
for i=1:3
    if strcmp(journey,'downlink')
        [posS,velS]=getSatPosVel(addSeconds(t,-dt),sids);
    else
        [posS,velS]=getSatPosVel(addSeconds(t,dt),sids);
    end
    dt=tof(posS,place);
end
end

function t=split2fields(toa)
t.d=day(toa,'dayofyear');
t.s=3600*toa.Hour+60*toa.Minute+toa.Second;
t.date=datetime(toa.Year,toa.Month,toa.Day);
end

function t1 = addSeconds(t,s)
t1=t;
t1.s=t.s+ s;
end

function [pos,vel]=getSatPosVel(toa,satIDs)
global jd2000;
global list;
UT1_UTC=-0.06;
TT_UTC=69.2;
noOfSats = length(satIDs);
pos = zeros(3,noOfSats);
vel = zeros(3,noOfSats);
for i=1:noOfSats
    satrec=list{satIDs(i)-400};
    epochDay=satrec.epochdays;
    ts = toa.s(i);
    d = toa.d(i);
    cdate = toa.date(i);
    tsince=(ts-(epochDay-floor(epochDay))*86400)/60 +(d-floor(epochDay))*1440;
    [~, pos1, vel1] = sgp4 (satrec,  tsince);
    jd=(ts/86400)+juliandate(cdate);
    jd_UT1 = jd + UT1_UTC/86400;
    jd_TT  = jd + TT_UTC/86400;
    ttt = (jd_TT-jd2000)/36525;
    [pos(:,i),vel(:,i),~]=teme2ecef(pos1',vel1',[0,0,0]',ttt,jd_UT1,0,0,0,2);
end
end

function fd=getDoppler(posS,velS,place,freq)
global LIGHTSPEED;
wvnum = freq/LIGHTSPEED;
dxyz=posS-place;
dr=sqrt(sum(dxyz.^2));
uvw=dxyz./dr;
vcomp=sum(uvw.*velS);
fd=-vcomp.*wvnum;
end

function dt=tof(pos1,pos2)
global LIGHTSPEED;
d=sqrt(sum((pos1-pos2).^2));
dt=d/LIGHTSPEED;
end

function [a]=getAngles(posS,place)
dxyz=posS-place;
r=sqrt(sum(dxyz.^2));
p=sqrt(sum(place.^2));
cosa=sum(dxyz.*place)./(r*p);
a=90-(acos(cosa)*180/pi);
end

function det=detDecesion(els)
pTable=[0,80,90,95,100,100,100,100,100,100,100,100,100,100,80,60,30,10]/100;
els(els<0)=0.1;
p=pTable(ceil(els/5));
x=rand(size(els));
det=x<p;
end