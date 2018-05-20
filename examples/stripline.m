
addpath(genpath([ pwd, '/..' ]));
%
% Compute per-length capacitance of a stripline.
%

% The geometry - trace above ground
mil2meter = 2.54e-5;
h = 10*mil2meter;     % height above ground
xt=12*mil2meter;    % trace width
yt=1*mil2meter;  % trace thickness
xg=12*mil2meter;    % width of the ground
yg=1*mil2meter;       % thickness of the ground
nxg=3;    % number of the segments along x in the ground
nyg=1;     % number of the segments along y in the ground
nxt=3;    % number of segments along x in the trace
nyt=1;    % number of segments along y in the trace

% trace
[ e, v ] = mkrect2d(xt, yt, nxt, nyt);
edges = e;
verts = v;

% lower ground
[ e, v ] = mkrect2d(xg, yg, nxg, nyg);
v = v + repmat([ 0 -(h+yg/2+yt/2) ], size(v, 1), 1);
edges = [ edges; e + size(verts, 1) ];
verts = [ verts; v ];

%% % upper ground
%% [ e, v ] = mkrect2d(xg, yg, nxg, nyg);
%% v = v + repmat([ 0  (h+yg/2+yt/2) ], size(v, 1), 1);
%% edges = [ edges; e + size(verts, 1) ];
%% verts = [ verts; v ];

% Find edges which belong to each of the conductors
c1 = 1:(nxt*2+nyt*2);
c2 = (nxt*2+nyt*2+1):size(edges, 1);
conductors = { c1 c2 };

%% plotmesh2d(edges, verts, { c2 }, 0)
%% xlim( [ -xg/2 xg/2 ] )
%% ylim( [ -xg/2 xg/2 ] )

nedges = size(edges, 1)
epsd = eps0;
epsout = repmat(epsd, nedges, 1);
epsin = 0*epsout;
C2 = extractc2(edges, verts, epsout, epsin, conductors)

%% C2 = [ 23.708e-12 -20.485e-12 ; -20.485e-12 23.708e-12 ]
%% L2 = [ 1851.926e-9    1600.152e-9 ; 1600.152e-9 1851.927e-9 ]

C = C2(1,1)%(C2(1,1)-C2(2,1))/2 % Mutual capacitance i.e. trace to ground
L = epsd*mu0/C
Z0=sqrt(L/C)
