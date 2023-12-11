function [TOT,FOA] = TRxOperation1(satIDs,ToA,FoT,TxSite,RxSite)

validChns=find(satIDs);
satIDs=satIDs(validChns);
noOfSats=length(validChns);
svID=satIDs-400; 

% t0=split2fields(ToA);
f_list=[1544.1e6,1544.9e6,1544.21e6];
%uplink
[posS,velS,dt1]=getPreciseSatPosVel_r(ToA+0.16/86400,svID,TxSite,'downlink');
% t1=addSeconds(t0,0.16);
% [posSx,velSx,dt1]=actualtof(t1,satIDs,TxSite,'uplink');
[posS1,velS1,~]=getPreciseSatPosVel_r(ToA+0.25/86400,svID,TxSite,'downlink');
% t2=addSeconds(t0,0.08);
% [posS1,velS1,~]=actualtof(t2,satIDs,TxSite,'uplink');%for doppler calculation
fd1=getDoppler(posS1,velS1,TxSite,FoT);
%downlink
dt2=tof(posS,RxSite);
cflag=floor(satIDs/100)-3;
freq_trns=f_list(cflag) - 406.05e6;
fc1=FoT+fd1+freq_trns;
fd2 = getDoppler(posS1,velS1,RxSite,fc1);
TOT=ToA-(dt1+dt2)/86400;
FOA=fc1+fd2;
end

function[posS,velS,dt]= actualtof(t,sids,place,journey)
dt=0;
for i=1:2
    if strcmp(journey,'downlink')
        [posS,velS]=getSatPosVel(addSeconds(t,-dt),sids);
    else
        [posS,velS]=getSatPosVel(addSeconds(t,dt),sids);
    end
    dt=tof(posS,place);
end
end

function[posS,velS,dt]=getPreciseSatPosVel_r(t,svID,place,journey)
global ephdata;
global eph_GLO;
noOfSat=length(t);
posS=zeros(3,noOfSat);
velS=zeros(3,noOfSat);
dt=zeros(1,noOfSat);
for i=1:noOfSat
    if svID(i)<100 %gallileo
        d=ephdata{svID(i)};
        [~,ind]=min(abs(seconds(t(i)-d.Time)));
        for j=1:2
            if strcmp(journey,'downlink')
                [posS_t,velS_t]=gnssconstellation(t(i)-dt(i)/86400,RINEXData=d(ind,:));
            else
                [posS_t,velS_t]=gnssconstellation(t(i)+dt(i)/86400,RINEXData=d(ind,:));
            end
            dt(i)=tof(posS_t.'*1e-3,place);
        end
        posS(:,i)=posS_t.'*1e-3;
        velS(:,i)=velS_t.'*1e-3;
    else %glonass
        %get the slot number
        tlist=[0,9,11,22];
        slno=tlist(svID(i)-100);
        t_jd=cal2jd_GT(t(i).Year,t(i).Month,t(i).Day+(t(i).Hour*3600+t(i).Minute*60+t(i).Second)/86400);
        t_jd = t_jd + 18/(60*60*24);
        [twk,Ttr,rollover]=jd2gps_GT(t_jd);
        [X,V]=SatPos_brdc_GLO(Ttr,slno+100,eph_GLO);
        posS(:,i)=X.'*1e-3;
        velS(:,i)=V.'*1e-3;
    end
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
