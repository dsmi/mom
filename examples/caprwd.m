
% Mutual capacity of two round wires in different dielectics,
% piecewise-linear approximation

addpath(genpath([ pwd, '/..' ]));

% The geometry - two round wires in two dielectrics
r = 5e-4;   % Radius of the wires
d = 1.5e-3; % Separation - distance between the centers
a = 1;      % Mesh accuracy multiplier
n = 4*a;    % Number of boundary elements around the wire
o = 0;      % dielectric boundary offset
l = 2e-2;   % dielectric boundary length
m = 20*a;   % Number of boundary elements in the dielectric boundary
epsd1 = eps0*1; % Upper dielectric
epsd2 = eps0*4; % Lower dielectric

[ e, v ] = mkcir2d(r, n);
v = v + repmat([ 0 -d/2 ], size(v, 1), 1);
edges = e;
verts = v;

[ e, v ] = mkcir2d(r, n);
v = v + repmat([ 0 d/2 ], size(v, 1), 1);
edges = [ edges; e + size(verts, 1) ];
verts = [ verts; v ];

[ e, v ] = mkline2d(l, m);
v = v + repmat([ o 0 ], size(v, 1), 1);
edges = [ edges; e + size(verts, 1) ];
verts = [ verts; v ];

% Find edges which belong to each of the conductors
conductor1 = find_edges2d( edges, verts, 0, -d/2, r*1.1 );
conductor2 = find_edges2d( edges, verts, 0, d/2, r*1.1 );
conductors = { conductor1 ; conductor2 };

% Edges forming the dielectric boundary
dieledges = find_eol2d( edges, verts, 0, 1, o );

nedges = size( edges, 1 );

% Outer and inner dielectric permittivities
epsout = zeros(nedges, 1);
epsout(conductor1) = epsd1;
epsout(conductor2) = epsd2;
epsout(dieledges) = epsd2;
epsin  = eps0 + 0*epsout;
epsin(dieledges) = epsd1;

bases = mkbases2d( edges );

%% plotmesh2d(edges,verts,{},1)
%% plotbases2d(edges,verts,bases)

C = extractc2( edges, verts, epsout, epsin, conductors );

Cmutual = (C(1,1)-C(2,1))/2 % Mutual capacitance

%% % piecewise-linear approximation
CL = extractc2l(edges, verts, mkbases2d( edges ), epsout, epsin, conductors );
Cmutlin = (CL(1,1)-CL(2,1))/2 % Mutual capacitance with pwl approximation

%% Cexact = 4.5526e-11; % calculated with detailed mesh
%% Cpwc_err = abs( Cexact - Cmutual ) / Cmutual
%% Cpwl_err = abs( Cexact - Cmutlin ) / Cmutual
