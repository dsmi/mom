%
% Compute resistance of a rectangular strip
%
addpath(genpath([ pwd, '/..' ]));

% The geometry - a strip
l = 5;
w = 1;
nl = 4*5;
nw = 4;
[ e, v ] = mkrect2d(l,w,nl,nw);

% Find ports
p1 = find_edges2d(e, v, -l/2, 0, w/2*(1+1e-8));
p2 = find_edges2d(e, v, l/2, 0, w/2*(1+1e-8));
ports = { p1' p2' };

plotmesh2d(e,v,ports,0);

Y = extracty2(e, v, ports, @intg_lapsl2d, @intg_lapdl2d)

R = Y(1,1)
R_test = w/l
