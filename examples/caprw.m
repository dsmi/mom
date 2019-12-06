
% Mutual capacity of two round/rectangular conductors,
% piecewise-linear approximation

addpath(genpath([ pwd, '/..' ]));

%% % The geometry - two round wires
%% r = 5e-4;   % Radius of the wires
%% d = 1.5e-3; % Separation - distance between the centers
%% n = 4;     % Number of boundary elements around the wire

%% [ e, v ] = mkcir2d(r, n);
%% v = v + repmat([ -d/2 0 ], size(v, 1), 1);
%% edges = e;
%% verts = v;

%% [ e, v ] = mkcir2d(r, n);
%% v = v + repmat([ d/2 0 ], size(v, 1), 1);
%% edges = [ edges; e + size(verts, 1) ];
%% verts = [ verts; v ];

% The geomentry -- two rectangular conductors
t = 1e-4;   % Thickness
w = 5e-4;   % Width
d = 1.5e-3; % Center-to-center separation
nt = 10; % Number of boundary elements along each edge
nw = 50; % Number of boundary elements along each edge

[ e, v ] = mkrect2d(w, t, nw, nt);
v = v + repmat([ -d/2 0 ], size(v, 1), 1);
edges = e;
verts = v;

[ e, v ] = mkrect2d(w, t, nw, nt);
v = v + repmat([ d/2 0 ], size(v, 1), 1);
edges = [ edges; e + size(verts, 1) ];
verts = [ verts; v ];

r = w; % So the bases lookup works correctly

%% plotmesh2d(edges,verts,{},0)
bases = mkbases2d( edges );

% Find edges which belong to each of the conductors
ce1 = find_edges2d( edges, verts, -d/2, 0, r*1.1 );
ce2 = find_edges2d( edges, verts, d/2, 0, r*1.1 );
cndedges = { ce1' ce2' };

% And bases which belong to each of the conductors
cb1 = find_bases2d( edges, verts, bases, -d/2, 0, r*1.1 );
cb2 = find_bases2d( edges, verts, bases,  d/2, 0, r*1.1 );
cndbases = { cb1' cb2' };

nedges = size( edges, 1 );

epsout = repmat(eps0, nedges, 1);
epsin = 0*epsout;
C = extractc2(edges, verts, epsout, epsin, cndedges);

Cmutual = (C(1,1)-C(2,1))/2 % Mutual capacitance

%% % Exact solution for round wires. d is the distance between centers.
%% % This assumes uniform current distribution, i.e. it works if
%% % the radius is much smaller than separation
%% c_per_l = pi*eps0/acosh(d/(2*r)) % Does not assume uniform current

%M = mkmommat2l( edges, verts, bases, @intg_lapsl2l );

% piecewise-linear approximation
CL = extractc2l(edges, verts, bases, eps0, cndbases);
Cmutlin = (CL(1,1)-CL(2,1))/2 % Mutual capacitance with pwl approximation

Cexact = 1.293119e-11
Cpwc_err = abs( Cexact - Cmutual ) / Cmutual
Cpwl_err = abs( Cexact - Cmutlin ) / Cmutual
