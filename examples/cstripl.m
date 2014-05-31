
addpath(genpath([ pwd, '/..' ]));
%
% Compute per-length capacitance of a microstrip.
%

% The geometry - trace above ground
h=5e-4;     % height above ground
xg=h*10;    % width of the ground
yg=h;       % thickness of the ground
nxg=200;    % number of the segments along x in the ground
nyg=20;     % number of the segments along y in the ground
xt=1e-2/32*4 % trace width
yt=1e-6;%xt/50;     % trace thickness
nxt=500;    % number of segments along x in the trace
nyt=1;       % number of segments along y in the trace

% trace
[ e, v ] = mkrect2d(xt, yt, nxt, nyt);
edges = e;
verts = v;

% another trace
%% [ e, v ] = mkrect2d(xt, yt, nxt, nyt);
%% v = v + repmat([ xt*2 0 ], size(v, 1), 1);
%% edges = [ edges; e + size(verts, 1) ];
%% verts = [ verts; v ];

% lower ground
[ e, v ] = mkrect2d(xg, yg, nxg, nyg);
v = v + repmat([ 0 -(h+yg/2+yt/2) ], size(v, 1), 1);
edges = [ edges; e + size(verts, 1) ];
verts = [ verts; v ];

% upper ground
[ e, v ] = mkrect2d(xg, yg, nxg, nyg);
v = v - repmat([ 0 -(h+yg/2+yt/2) ], size(v, 1), 1);
edges = [ edges; e + size(verts, 1) ];
verts = [ verts; v ];

% Find edges which belong to each of the conductors
c1 = 1:(nxt*2+nyt*2);
c2 = (nxt*2+nyt*2+1):size(edges, 1);
conductors = { c1 c2 };

%plotmesh2d(edges, verts, { c2 }, 1)

nedges = size(edges, 1);

epsout = repmat(eps0, nedges, 1);
epsin = 0*epsout;
C2 = extractc2(edges, verts, epsout, epsin, conductors);
C=C2(1,1)
L=eps0*mu0/C
