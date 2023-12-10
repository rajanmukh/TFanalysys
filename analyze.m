 global eph_GLO;
if ~exist('initialized','var')
    addpath([pwd,'\sgp4']);
    initializeRecord();
    readrinex('IISC00IND_R_20233380000_01D_EN.rnx');
    eph_GLO=read_nav_glo('IISC00IND_R_20233380000_01D_RN.rnx',18);
    initialized = true;
end
fileID=fopen('beacondata_LE_2023_12_04.txt');
load('prns.mat')
present_hour = 0;
prev_hour = 0;
present_day=0;
prev_day=0;
f_list=[1544.1e6,1544.9e6,1544.21e6];
RxSite=lla2ecef([13.036,77.5124,930])'*1e-3;%Bangalore
TxSite=lla2ecef([43.56,1.48,214])'*1e-3;%France
FoT=406.022e6;
REFID = '9C6000000000001';
noB=2880;
TOA_m=repmat(NaT,7,noB);
terr=zeros(7,noB);
ferr=zeros(7,noB);
sats = zeros(7,noB);
idx=1;

while(1)
    ca=fgetl(fileID);
    if isempty(ca)        
        continue;
    end
    if ~ischar(ca)
        break;%end of file
    end
    
    str=string(ca);
    ss=split(ca,',');
    antmarks=str2double(ss(13:19));
    snames=ss([25,46,67,88,109,130,151]);
    ants=find(antmarks==1 & ~strcmp(snames,'DEFAULT'));
    nos = length(ants);
    msgs = cell(1,nos);
    brates= cell(1,nos);
    toa_m = repmat(NaT,1,nos);
    tmstamps = cell(1,nos);
    foas = zeros(1,nos);
    CNRs = zeros(1,nos);
    SIDs = zeros(1,nos);
    pdf1errs=cell(1,nos);
    pdf2errs=cell(1,nos);
    for i=1:nos
        chn=ants(i);
        ii=23+(chn-1)*21; 
        
        sname = char(ss(ii+2));
        SIDs(i)=getSID(sname,prns);
        tmstamps{i} = char(ss(ii+3));
        brates{i} = char(ss(ii+4));
        id=char(ss(ii+7));
        otherID=~strcmp(id,REFID);
        if otherID
            break;
        end
        msg=char(ss(ii+8));
        
        msgs{i}=msg;
        CNRs(i)=str2double(ss(ii+9));
        foas(i)=str2double(ss(ii+10));
        pdf1errs{i}=char(ss(ii+11));
        pdf2errs{i}=char(ss(ii+12));
        toa = ss(ii+13:ii+20);
        if isempty(toa(8))
            continue;
        end
        TOA=[num2str(str2double(toa(1))+2000),'-',num2str(str2double(toa(2)),'%03d'),' ',num2str(str2double(toa(3)),'%02d'),':',num2str(str2double(toa(4)),'%02d'),':',num2str(str2double(toa(5)),'%02d'),':',num2str(str2double(toa(6)),'%03d'),num2str(str2double(toa(7)),'%03d'),num2str(str2double(toa(8)),'%03d')];
        toa_m(i) = datetime(TOA,'InputFormat','uuuu-DDD HH:mm:ss:SSSSSSSSS'); 
        present_hour = str2double(toa(3));
        present_day = str2double(toa(2));
        if present_hour ~= prev_hour
            readtle(toa_m(i));
        end           
    end
    if otherID
        continue;
    end
    if nos>0
        %find out max CNR packet
        [~,ind]=max(CNRs);
        data = hexToBinaryVector(msgs{ind},144);
        jday=100*binaryVectorToDecimal(data(113:116))+10*binaryVectorToDecimal(data(117:120))+binaryVectorToDecimal(data(121:124));
        hr=10*binaryVectorToDecimal(data(125:128))+binaryVectorToDecimal(data(129:132));
        mn=10*binaryVectorToDecimal(data(133:136))+binaryVectorToDecimal(data(137:140));   
        sec=10*binaryVectorToDecimal(data(141:144));
        ToT=datetime(2023,12,4,hr,mn,sec);
        cflag=floor(SIDs/100)-3;
%         if ~any(cflag==2)
%             continue;
%         end
        freq_trns=f_list(cflag) - 406.05e6;
        foa_m = foas+53.1311e3+f_list(cflag)-1e5;        
        [toa_e,foa_e] = TRxOperation(SIDs,ToT,FoT,TxSite,RxSite); 
        terr(ants,idx)=seconds(toa_m-toa_e);
        ferr(ants,idx)=foa_m-foa_e;
        TOA_m(ants,idx)=toa_m;
        sats(ants,idx)=SIDs;
        idx=idx+1;idx
    end    
    prev_hour = present_hour;
    prev_day = present_day;
end
fclose(fileID);
chn=1;
valid=~isnat(TOA_m);
ts=TOA_m(chn,valid(chn,:));
tr=terr(chn,valid(chn,:));
fr=ferr(chn,valid(chn,:));
sr=sats(chn,valid(chn,:));
ix=find(abs(fr+13)>5);
ts(ix)=[];
tr(ix)=[];
fr(ix)=[];
sr(ix)=[];




