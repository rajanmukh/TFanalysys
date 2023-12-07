function fd=getDoppler(posS,velS,place,freq)
LIGHTSPEED = 299792.458;%km/sec
dxyz=posS-place;
dr=sqrt(sum(dxyz.^2));
uvw=dxyz./dr;
vcomp=sum(uvw.*velS);
fd=-vcomp.*freq/LIGHTSPEED;
end

