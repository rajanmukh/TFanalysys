function [loc,err,antsV,sInfo] = computeLocation(toas,foas,cnrs,satIDs)
loc=[];
err=[];
sInfo=[];
%constant
%station location
BANGALORE=1e3*[1.344164600515364   6.068648167092199   1.429495327500622]';
% STATIONARY=[0,0,0]';
f_list=[1544.1e6,1544.9e6,1544.21e6];
noOfChns=length(satIDs);
cflag=floor(satIDs/100)-3;
freq_trns=f_list(cflag) - 406.05e6;
fc = foas+53.1311e3+f_list(cflag)-1e5;
t0=split2fields(toas);
% [toa,d,date,fc,freq_trns]=readSatParamsSGP(datetoa,foa,satIDs);
[stdtoa,stdfoa]=estimateMeasError(cnrs);
t1=addSeconds(t0,0.16);
[posS,velS,dt]=actualtof(t1,satIDs,BANGALORE,'downlink');
sInfo.satPos=posS;
sInfo.satVel=velS;
sInfo.toff=dt;
sInfo.upTOA = toas - dt/86400;
t=t0.s-dt;%onboard transmit/receive time
t2=addSeconds(t0,0.08);
[posS1,velS1,~]=actualtof(t2,satIDs,BANGALORE,'downlink');%for doppler calculation
fd=getDoppler(posS1,velS1,BANGALORE,fc);
sInfo.foff=fd;
fc1=fc-fd-freq_trns;
sInfo.upFOA=fc1;
noOfSats=length(satIDs);
if noOfSats>=2    
    G=firstGuess(posS,t);
    G(5)=mean(fc1);
    for i=1:15
        %     [posS1,dt2]=actualtof2(posS,G(1:3));
        [F,D]=FDcreatorTD(posS,t,posS1,velS1,fc1,G,stdtoa,stdfoa,noOfSats);
        %     [F,D]=FDcreator3(posS1,velS1,fc1,G);
        %     [F,D]=FDcreator1(posS,t,G);
        del=D\F;
%         [del,~,resd]=lscov(D,F);
        del(4)=del(4)*1e-3;
        G=G-del;
    end
    resd = norm(F);

    freqOutOfBounds = G(5)<406.01e6 || G(5)>406.09e6;    

    if resd>100 || freqOutOfBounds       
        errorDetected=true;        
        %try to identify wrong channel
%         [~,tryOrder] = sort(abs(F(1:noOfChns)),'descend');
        [G,D,antsV,errorEliminated] = tryElimination(satIDs,posS,t,posS1,velS1,fc1,stdtoa,stdfoa);        
    else        
        errorDetected=false;
        antsV=ones(1,noOfChns,'logical');
    end

    if ~errorDetected || (errorDetected && errorEliminated)
        location_est=G(1:3);
        % t(:)=G(4);
        [posS2,velS2,~]=actualtof(t2,satIDs,location_est,'uplink');
        fd1 = getDoppler(posS2,velS2,location_est,fc1);
        ft=fc1-fd1;
        llaPos=ecef2lla(1e3*G(1:3)');
        loc.lat=llaPos(1);
        loc.lon=llaPos(2);
        loc.alt=llaPos(3)/1e3;
        [jdop,elp]=computeDOP(D,loc.lat*pi/180,loc.lon*pi/180);
        err.EHE = 2.5*jdop;
        
        err.ellipse=elp;
        if errorDetected
            stdtoa = stdtoa(antsV);
            ft = ft(antsV);
            noOfSats = sum(antsV);
        end
        sigt=median(stdtoa);
        sInfo.jdop=jdop/sigt;
        sInfo.ft=mean(ft);
        if noOfSats == 3
            sInfo.solMethodology = 'P2Dd-Single';
        else
            sInfo.solMethodology = 'P3D-Single';
        end        
    end
else
    antsV=ones(1,noOfChns,'logical');%consistency cannot be checked for lack of sats
end
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

function [pos]=adjustRotation(xyz,dt)
omega_dot_earth = 7.2921151467e-5; %(rad/sec)
ths=omega_dot_earth*dt;
pos=zeros(size(xyz));
for i=1:length(dt)
    th=ths(i);
    R=[cos(th) sin(th) 0; -sin(th) cos(th) 0;0 0 1];
    pos(:,i)=R*xyz(:,i);
end
end

function dt=tof(pos1,pos2)
global LIGHTSPEED;
d=sqrt(sum((pos1-pos2).^2));
dt=d/LIGHTSPEED;
end

function [F,D] = FDcreator(posS,t,G)
global LIGHTSPEED;
xyz=G(1:3);
tg=G(4);
R=sqrt(sum((xyz-posS).^2));
F=(R-LIGHTSPEED*(t-tg))';
D=zeros(length(t),4);
D(:,1:3)=((1./R).*(xyz-posS))';
D(:,4)=1e-3*LIGHTSPEED;
end

function [F,D] = FDcreator1(posS,t,G)
global LIGHTSPEED;
xyz=G(1:3);
tg=G(4);
R=sqrt(sum((xyz-posS).^2));
R(end+1)=sqrt(sum(xyz.^2));
f=1/298.257223560;
a=6378.137;
sinth2=(xyz(3))^2/(R(end)^2);
lrad=a*(1-f*sinth2);
obs_range=[LIGHTSPEED*(t-tg) lrad];
F=(R-obs_range)';
D=zeros(length(t)+1,4);
D(1:length(t),1:3)=((1./R(1:length(t))).*(xyz-posS))';
D(1+length(t),1:3)=((1./R(1+length(t))).*xyz)';
D(1:length(t),4)=1e-3*LIGHTSPEED;
end

function [F,D] = FDcreatorTD(posS,t,posS1,velS1,freq,G,stdtoa,stdfoa,noOfSats)
if noOfSats == 3 || noOfSats == 2
    [F1,D1]=FDcreator1(posS,t,G(1:4));
    stdtoa=[stdtoa;0.5];
else
    [F1,D1]=FDcreator(posS,t,G(1:4));
end
[F2,D2]=FDcreator3(posS1,velS1,freq,G([1,2,3,5]));
F=[F1./stdtoa;F2./stdfoa];
D=zeros(length(F),5);
D(1:length(F1),1:4)=D1./stdtoa;
D(length(F1)+1:length(F),[1:3,5])=D2./stdfoa;
end

function [F,D] = FDcreatorT(posS,t,posS1,velS1,freq,G,stdtoa,stdfoa,noOfSats)
if noOfSats == 3 || noOfSats == 2
    [F1,D1]=FDcreator1(posS,t,G(1:4));
    stdtoa=[stdtoa;0.5];
else
    [F1,D1]=FDcreator(posS,t,G(1:4));
end
F=F1./stdtoa;
D=D1./stdtoa;
end

function [F,D] = FDcreator3(posS,velS,freq,G)
global LIGHTSPEED;
noOfSat=length(freq);
F=zeros(noOfSat,1);
xyz=G(1:3);
fg=G(4);
dxyz=posS-xyz;
dr=sqrt(sum(dxyz.^2));
uvw=dxyz./dr;
vcomp=sum(uvw.*velS);
wvlen = LIGHTSPEED./freq;
fd=freq-fg;
F(1:noOfSat)=(vcomp+wvlen.*fd)';

D=zeros(noOfSat,4);
D(:,1:3) = ((1./dr).*(-velS+uvw.*vcomp))';
D(:,4) = -wvlen';
end

function [G] = firstGuess(pos,t)
%FIRSTGUESS Summary of this function goes here
%   Detailed explanation goes here
EARTHCENTER=[1000;6000;1000];
dt=tof(pos,EARTHCENTER);
G=[EARTHCENTER;mean(t-dt)];
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

function [posS,dt]=actualtof2(posS,place)
for j=1:3
    dt=tof(posS,place);
    [posS] = adjustRotation(posS,-dt);
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

function [JDOP,elpse] = computeDOP(D,lat,lon)
R=[-sin(lon) cos(lon) 0;-sin(lat)*cos(lon) -sin(lat)*sin(lon) cos(lat);cos(lat)*cos(lon) cos(lat)*sin(lon) sin(lat)];
P = inv(D'*D);
Q=R*P(1:3,1:3)*R';
JDOP=sqrt(Q(1,1)+Q(2,2));
[V,D]=eig(Q(1:2,1:2));
[ab,ind]=sort(sqrt(diag(D)),'descend');
maxind=ind(1);
A=atan(V(2,maxind)/V(1,maxind))*180/pi;
if A<0
    A=360+A;
end
elpse=[1.4*ab(1),1.4*ab(2),A];
end

function [terrstd,ferrstd] = estimateMeasError(cbn0)
terrstd=0.3*(15*2.^((-cbn0+35)/6))';
ferrstd=0.7e-3*(0.2*2.^((-cbn0+35)/6)+1)';
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

