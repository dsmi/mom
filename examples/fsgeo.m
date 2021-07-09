
% Load fieldsolver-saved mesh and solve it

addpath(genpath([ pwd, '/..' ]));

fname = 'd:/fieldsolver_geo_7.txt';
[ verts, edges, bases, erout, erin, conductors ] = load_mesh( fname );

% Flip y -- this fixes the normals direction and makes the geometry upside up
verts = [ verts(:,1) -verts(:,2) ];

% conductor-to-dielectric edges
cndedges = cell2mat(conductors);

% dielectric-to-dielectric edges
dieledges = find( ~ismember( 1:size(edges,1), cndedges ) );

% Bases forming the dielectric boundaries
dielectric_bases = bases( find( ~ismember( bases(:,1), cndedges ) ), : );

%% plotbases2d(edges,verts,dielectric_bases)

% Scale dielectric verts -- temporary experiment
dieledges = edges( find( ~ismember( transpose( 1:size( edges, 1 ) ), cndedges ) ), : );
verts( [ dieledges( :, 1 ) ; dieledges( :, 2 ) ], : ) = verts( [ dieledges( :, 1 ) ; dieledges( :, 2 ) ], : ) * 1.0001;

%% plotmesh2d(edges,verts,{ },0)

% Omit the losses
epsout = eps0 * real( erout );
epsin = eps0 * real( erin );

C3 = extractc2( edges, verts, epsout, epsin, conductors );
C = cperlen( C3 )

% piecewise-linear approximation
Cl3 = extractc2l( edges, verts, bases, epsout, epsin, conductors );
Cl = cperlen( Cl3 )

%% Gap test, 
%%   C = 6.22862e-11 (fieldsolver with acc = 8)
%%   C = 6.2224e-11 (fieldsolver with acc = 2)
%%  Accuracy -6, 2243 unknowns
%%   C  = 6.225673372246184e-11
%%   Cl = 6.227175075429566e-11
%%  Accuracy -3, 274 unknowns
%%   C  = 6.196087085608778e-11
%%   Cl = 6.218993201307802e-11


%% Symmetric test
%% C  =  7.24006e-11 (accuracy = 8, 14807 unknowns)
%% Cl =  0.000000000072467
%% Cl =  0.000000000072467 (7397 unknowns)

%% Result for 'mytest.inp' with accuracy = 8
%% Capacitance matrix [pF/m]:
%%         1        2
%%  1   74.0380  -9.8211
%%  2   -9.8211  74.0404

%% %% Result for FS_Error50.inp, accuracy = 7, 24890 unknowns
%% Cx = ...    
%%    [ 89.7083  -5.8118  -0.4911  -0.1018 ;...
%%      -5.8118  89.7736  -1.8050  -0.1548 ;...
%%      -0.4911  -1.8050  99.7459  -0.3140 ;...
%%      -0.1018  -0.1548  -0.3140  99.6522 ]*1e-12;
%% %% Accuracy = 2, 769 unknowns
%% Cfs = ...
%%   [ 89.5015  -5.8671  -0.4865  -0.0818 ;...
%%    -5.8671  89.5671  -1.8138  -0.1215 ;...
%%    -0.4865  -1.8138  99.5340  -0.3243 ;...
%%    -0.0818  -0.1215  -0.3243  99.4512 ]*1e-12;

%% Cl_err = norm(Cl-Cx)
%% Cfs_err = norm(Cfs-Cx)

%% 'exact' result with very fine mesh
%%        1        2
%% 1   78.1523 -19.7677
%% 2  -19.7677  78.1601
