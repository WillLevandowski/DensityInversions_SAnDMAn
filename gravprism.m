%---------------------------------------------------------------------
function dg = gravprism(drho,dx1,dx2,dy1,dy2,dz1,dz2)
%
% gravitational attraction due to "m" prisms at "n" observation point
% x1,x2,y1,y2,z1,z2 are coordinates of edges of prisms relative to 
% observation points.  They are m x n matrices
%---------------------------------------------------------------------

gam=(6.6732e-11)*1e5;             % mGal m^2/kg
g=0;

R111=sqrt(dx1.^2+dy1.^2+dz1.^2); 
R112=sqrt(dx2.^2+dy1.^2+dz1.^2);
R121=sqrt(dx1.^2+dy2.^2+dz1.^2);
R122=sqrt(dx2.^2+dy2.^2+dz1.^2);
R211=sqrt(dx1.^2+dy1.^2+dz2.^2); 
R212=sqrt(dx2.^2+dy1.^2+dz2.^2);
R221=sqrt(dx1.^2+dy2.^2+dz2.^2);
R222=sqrt(dx2.^2+dy2.^2+dz2.^2);

g111=-[dz1.*atan((dx1.*dy1)./(dz1.*R111))-dx1.*log(R111+dy1)-dy1.*log(R111+dx1)];
g112=+[dz1.*atan((dx2.*dy1)./(dz1.*R112))-dx2.*log(R112+dy1)-dy1.*log(R112+dx2)];
g121=+[dz1.*atan((dx1.*dy2)./(dz1.*R121))-dx1.*log(R121+dy2)-dy2.*log(R121+dx1)];
g122=-[dz1.*atan((dx2.*dy2)./(dz1.*R122))-dx2.*log(R122+dy2)-dy2.*log(R122+dx2)];

g211=+[dz2.*atan((dx1.*dy1)./(dz2.*R211))-dx1.*log(R211+dy1)-dy1.*log(R211+dx1)];
g212=-[dz2.*atan((dx2.*dy1)./(dz2.*R212))-dx2.*log(R212+dy1)-dy1.*log(R212+dx2)];
g221=-[dz2.*atan((dx1.*dy2)./(dz2.*R221))-dx1.*log(R221+dy2)-dy2.*log(R221+dx1)];
g222=+[dz2.*atan((dx2.*dy2)./(dz2.*R222))-dx2.*log(R222+dy2)-dy2.*log(R222+dx2)];

dg=drho.*gam.*(g111+g112+g121+g122+g211+g212+g221+g222);