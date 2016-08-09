
addpath(genpath([ pwd, '/..' ]));
%
% Compute per-length capacitance of a stripline.
%

% The geometry - trace above ground
h=3e-5;     % height above ground
xg=h*50;    % width of the ground
yg=h;       % thickness of the ground
nxg=500;    % number of the segments along x in the ground
nyg=100;     % number of the segments along y in the ground
xt=5e-5;    % trace width
yt=1.5e-5;  % trace thickness
nxt=300;    % number of segments along x in the trace
nyt=50;    % number of segments along y in the trace

% trace
[ e, v ] = mkrect2d(xt, yt, nxt, nyt);
edges = e;
verts = v;

% lower ground
[ e, v ] = mkrect2d(xg, yg, nxg, nyg);
v = v + repmat([ 0 -(h+yg/2+yt/2) ], size(v, 1), 1);
edges = [ edges; e + size(verts, 1) ];
verts = [ verts; v ];

% Find edges which belong to each of the conductors
c1 = 1:(nxt*2+nyt*2);
c2 = (nxt*2+nyt*2+1):size(edges, 1);
conductors = { c1 c2 };

%plotmesh2d(edges, verts, { c2 }, 0)

nedges = size(edges, 1);
epsout = repmat(eps0, nedges, 1);
epsin = 0*epsout;
C2 = extractc2(edges, verts, epsout, epsin, conductors);

C = (C2(1,1)-C2(2,1))/2 % Mutual capacitance i.e. trace to ground
L=eps0*mu0/C
Z0=sqrt(L/C)
