function readtle(dtime)
opsmode='a';
typerun='m';
typeinput='m';
whichconst=84;
global list;
list=cell(1,300);
load('prns.mat','prns');
files=dir('TLE\tle*');
if ~isempty(files)
%     infile=-1;
    filenames = arrayfun(@(x) string(x.name(5:end)),files);
    dtarr=datetime(filenames,'InputFormat','DDD_uuuu');
    [~,ind]=min(abs(seconds(dtime-dtarr)));
    infilename=['TLE\',files(ind).name];
    infile = fopen(infilename, 'r');    
    if infile > 0
        while true
            longstr=fgets(infile);
            if length(longstr)<2
                break;
            else
                if startsWith(longstr,'BEIDOU')
                    continue;
                elseif startsWith(longstr,'COSMOS')
                    satno=str2double(extractBetween(longstr,10,11));
                elseif startsWith(longstr,'GSAT')
                    satno=str2double(extractBetween(longstr,6,8));
                else
                    continue;
                end
                if satno>0 && satno<=300
                    satID=prns(satno);
                end
                if satID > 400
                    longstr1=fgets(infile);
                    longstr2=fgets(infile);
                    % sgp4fix additional parameters to store from the TLE
                    rec.classification = 'U';
                    rec.intldesg = '        ';
                    rec.ephtype = 0;
                    rec.elnum   = 0;
                    rec.revnum  = 0;
                    
                    [~,~,~, rec] = twoline2rv( ...
                        longstr1, longstr2, typerun, typeinput, opsmode, whichconst);
                    
                    list{satID-400}=rec;
                end
            end
            
        end
    end
    fclose(infile);    
end
end


