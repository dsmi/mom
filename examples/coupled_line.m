
addpath(genpath([ pwd, '/..' ]));
%
% Compute per-length capacitance of a coupled two-conductor line.
%

% The geometry - trace above another trace
mil2meter = 2.54e-5;
h  = 10*mil2meter    % vertical separation
xt = 12*mil2meter    % trace width
yt = 1.35*mil2meter  % trace thickness
nxg=1;    % number of the segments along x in the ground
nyg=1;    % number of the segments along y in the ground
nxt=50;    % number of segments along x in the trace
nyt=7;    % number of segments along y in the trace

% trace 1
[ e1, v1 ] = mkrect2d(xt, yt, nxt, nyt);
v1 = v1 + repmat([ 0 -(h+yt)/2 ], size(v1, 1), 1);

% trace 2
[ e2, v2 ] = mkrect2d(xt, yt, nxt, nyt);
v2 = v2 + repmat([ 0 (h+yt)/2 ], size(v2, 1), 1);

% join meshes
edges = [ e1 ; e2 + size(v1, 1) ];
verts = [ v1 ; v2 ];

% Find edges which belong to each of the conductors
c1 = 1:size(e1, 1);
c2 = (size(e1, 1)+1):(size(e1, 1)+size(e2, 1));
conductors = { c1 c2 };

plotmesh2d(edges, verts, { c2 }, 0)
xlim( [ -2*h 2*h ] )
ylim( [ -2*h 2*h ] )

nedges = size(edges, 1)
epsd = 4.3*eps0;
epsout = repmat(epsd, nedges, 1);
epsin = 0*epsout;
C2 = extractc2(edges, verts, epsout, epsin, conductors)

C = cperlen(C2)

%% L = epsd*mu0/C
%% Z0=sqrt(L/C)
