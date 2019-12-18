function [ verts, edges, bases, epsout, epsin, conductors ] = load_mesh( fname )
% [ verts, edges, bases, epsout, epsin, conductors ] = load_mesh( fname )
%
% Read mesh for extractc2l from a text file. Sample test file can be
% found in the tests area.
% No error checking is done.
%

f = fopen( fname, 'rt' );

% Vertices
nv = fscanf( f, '%i', 1 );
verts = transpose( fscanf( f, '%f', [ 2, nv ] ) );

% Edges
ne = fscanf( f, '%i', 1 );
edges = transpose( fscanf( f, '%i', [ 2, ne ] ) );

% Bases
nb = fscanf( f, '%i', 1 );
bases = transpose( fscanf( f, '%i', [ 2, nb ] ) );

% Permittivity outside
neo = fscanf( f, '%i', 1 );
epsoutri = transpose( fscanf( f, '%f %fi', [ 2, neo ] ) );
epsout = epsoutri(:,1) + i*epsoutri(:,2);

% Permittivity inside
nei = fscanf( f, '%i', 1 );
epsinri = transpose( fscanf( f, '%f %fi', [ 2, nei ] ) );
epsin = epsinri(:,1) + i*epsinri(:,2);

% Conductors
conductors = {};
nc = fscanf( f, '%i', 1 );
for cidx = 1:nc
    nce = fscanf( f, '%i', 1 );
    ce = fscanf( f, '%i', nce );
    conductors{ cidx, 1 } = ce;
end

fclose( f );

