function [a,b,r]=getAngles(posS,place)
dxyz=posS-place;
r=sqrt(sum(dxyz.^2));
p=sqrt(sum(place.^2));
s=sqrt(sum(posS.^2));
cosa=sum(dxyz.*place)./(r*p);
a=90-(acos(cosa)*180/pi);
cosb=sum(dxyz.*posS)./(r.*s);
b=acos(cosb)*180/pi;
end
