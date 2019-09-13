function [ pN,tN ] = intpol( p,t,N )
%Interpolates between points
% You get to pick how many points, given by input N
% Dr. Muhammar Qureshi 12/15/2016

t0 = t(1);
tend = t(end); 
m = N-1;

n = length(p)-1;
deltat = (tend - t0)/n;
told = t0:deltat:tend;

deltatnew = (tend - t0)/m;
tnew = t0:deltatnew:tend;

pN = interp1(t,p,tnew,'pchip');
% qN = interp1(told,q,tnew,'spline');
tN = interp1(told,t,tnew,'linear');

end

