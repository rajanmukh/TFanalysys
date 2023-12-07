function dt=tof(posS,place)
LIGHTSPEED = 299792.458;%km/sec
dxyz=posS-place;
r=sqrt(sum(dxyz.^2));
dt=r/LIGHTSPEED;
end
