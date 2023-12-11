function SID=getSID(snm,prnlist)
if snm(1)=='G'
    idx=str2double(extractAfter(snm,5));
elseif snm(1)=='C'
    idx=str2double(extractAfter(snm,9));
else
    idx=2;
end
SID=prnlist(idx);
end


