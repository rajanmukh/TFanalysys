function satname = getSatName(sid,list)
satname='';
d=find(list==sid);
if sid>500
    satname = ['COSMOS 25',num2str(d,'%02d')];
elseif sid>400
    satname = ['GSAT0',num2str(d)];
end
end
