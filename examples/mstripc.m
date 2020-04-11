
%
% Compute per-length capacitance of a microstrip line with 
% a dielectric layer between the line and the ground
%
addpath(genpath([ pwd, '/..' ]));

% The geometry - microstrip above ground, dielectric in between
lu = 2.54e-5; % length unit, mil-to-meter
w  = 5.0*lu;  % microstrip width
t  = 1.35*lu; % microstrip thickness
a  = w*40;    % dielectric and ground width
b  = 7.0*lu;  % dielectric thickness
d  = (a/2-w/2);  % aux, with of a half of the dielectric interface
eps2 = eps0*4.2; % Dielectric
n  = 1;       % Mesh detail level -- edges along the short trace edge
nw = ceil(n*w/t); % Mesh parameters below 
nt = n;           % 
na = nw*4; %
nd = nw*2; %

% Ground!
[ e, v ] = mkline2d( a, na );
v = v + repmat( [ 0 -b ], size(v, 1), 1 );
edges = e;
verts = v;

% Trace
[ e, v ] = mkrect2d( w, t, nw, nt );
v = v + repmat( [ 0 t/2 ], size(v, 1), 1 );
[ edges, verts ] = joinmeshes2d( edges, verts, e, v );

% Find edges which belong to each of the conductors
c0 = find_eol2d( edges, verts, 0, 1, -b );
c1 = find_edges2d( edges, verts, 0, 0, w/2+t/2 );
cnd_edges = { c0 ; c1 };

% Now add the dielectric interface
[ e, v ] = mkrevlogline2d( d, nd );
v = v + repmat( [ -d/2-w/2, 0 ], size(v, 1), 1 );
[ edges, verts ] = joinmeshes2d( edges, verts, e, v );

[ e, v ] = mklogline2d( d, nd );
v = v + repmat( [ d/2+w/2, 0 ], size(v, 1), 1 );
[ edges, verts ] = joinmeshes2d( edges, verts, e, v );

%% plotmesh2d( edges, verts, cnd_edges', 1 );

% Dielecric permeabilities in and out
epsout = eps0 * ones( size( edges, 1 ), 1 );
epsin = eps0 * ones( size( edges, 1 ), 1 );

% Dielectric-to-air interface, and the bottom surface of the trace
d2a = find_eol2d( edges, verts, 0, 1, 0 );
epsout( d2a ) = eps2;

% Ground plane
gpl = find_eol2d( edges, verts, 0, 1, -b );
epsin( gpl ) = eps2;

bases = mkbases2d( edges );

% hack -- add base joining the dielectric surfaces
bases = [ bases ; [ size( edges, 1 ) - nd, size( edges, 1 ) - nd + 1 ] ];

% All the edges of the conductor surfaces
allcndedges = cell2mat( cnd_edges );

%% % Bases forming the dielectric boundaries
%% dielectric_bases = bases( find( ~ismember( bases(:,1), allcndedges ) ), : );

%% plotbases2d(edges,verts,dielectric_bases)


C = cperlen( extractc2( edges, verts, epsout, epsin, cnd_edges ) )
Cl = cperlen( extractc2l( edges, verts, bases, epsout, epsin, cnd_edges ) )
