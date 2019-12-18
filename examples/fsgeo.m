
% Load fieldsolver-saved mesh and solve it

addpath(genpath([ pwd, '/..' ]));

fname = 'fieldsolver_geo_4.txt';
[ verts, edges, bases, erout, erin, conductors ] = load_mesh( fname );

% Flip y -- this fixes the normals direction and makes the geometry upside up
verts = [ verts(:,1) -verts(:,2) ];

% conductor-to-dielectric edges
cndedges = cell2mat(conductors);

% dielectric-to-dielectric edges
dieledges = find( ~ismember( 1:size(edges,1), cndedges ) );

%% % Bases forming the dielectric boundaries
%% dielectric_bases = bases( find( ~ismember( bases(:,1), cndedges ) ), : );

%% plotbases2d(edges,verts,dielectric_bases)

%% plotmesh2d(edges,verts,{},1)

% Omit the losses
epsout = eps0 * real( erout );
epsin = eps0 * real( erin );

C3 = extractc2( edges, verts, epsout, epsin, conductors );
C = cperlen( C3 )

% piecewise-linear approximation
Cl3 = extractc2l( edges, verts, bases, epsout, epsin, conductors );
Cl = cperlen( Cl3 )

 %% 'exact' result with very fine mesh
 %%        1        2
 %% 1   78.1523 -19.7677
 %% 2  -19.7677  78.1601

