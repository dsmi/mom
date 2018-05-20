
addpath(genpath([ pwd, '/..' ]));

%
% Compute per-length capacitance of a coupled two-conductor line.
%

% The geometry - trace above another trace
mil2meter = 2.54e-5;
h  = 10*mil2meter    % vertical separation
w  = 12*mil2meter    % trace width
gw = 200*mil2meter    % ground width
t = 1.35*mil2meter  % trace thickness
nw=50;    % number of segments along x in the trace
nt=7;     % number of segments along y in the trace
ng=1000;  % number of segments along y in the ground

% trace 1
[ e1, v1 ] = mkrect2d(w, t, nw, nt);
v1 = v1 + repmat([ 0 (h+t) ], size(v1, 1), 1);

% trace 2
[ e2, v2 ] = mkrect2d(w, t, nw, nt);

% trace 3
[ e3, v3 ] = mkrect2d(gw, t, ng, nt);
v3 = v3 + repmat([ 0 -(h+t) ], size(v3, 1), 1);

% join meshes
edges = [ e1 ; e2 + size(v1, 1) ; e3 + size(v1, 1) + size(v2, 1) ];
verts = [ v1 ; v2 ; v3 ];

% Find edges which belong to each of the conductors
c1 = 1:size(e1, 1);
c2 = (size(e1, 1)+1):(size(e1, 1)+size(e2, 1));
c3 = (size(e1, 1)+size(e2, 1)+1):size(edges, 1);
conductors = { c3 c1 c2 }; % first is ground

%% plotmesh2d(edges, verts, { conductors{1} }, 0)
%% xlim( [ -gw/2 gw/2 ] )
%% ylim( [ -gw/2 gw/2 ] )

nedges = size(edges, 1)
epsd = 4.3*eps0;
epsout = repmat(epsd, nedges, 1);
epsin = 0*epsout;
C2 = extractc2(edges, verts, epsout, epsin, conductors)

C = cperlen(C2)

L = epsd*mu0*inv(C)
%% Z0=sqrt(L/C)
