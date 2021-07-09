function test_find_bases2d
% Test find_edges2d function, which finds bases with both of the support
% edges in the given circle
%

% The testing geometry - a pair of circles
r = 5e-4;   % Radius
d = 1.5e-3; % Separation
n = 6;      % Number of the segments in each of the circles

[ e, v ] = mkcir2d(r, n);
v = v + repmat([ -d/2 0 ], size(v, 1), 1);
edges = e;
verts = v;

[ e, v ] = mkcir2d(r, n);
v = v + repmat([ d/2 0 ], size(v, 1), 1);
edges = [ edges; e + size(verts, 1) ];
verts = [ verts; v ];

%% plotmesh2d(edges,verts,{},0)
bases = mkbases2d( edges );

% This is expected to find bases belonging to the first circle
cb1 = find_bases2d( edges, verts, bases, -d/2, 0, r*1.1 );
assertEquals( (1:6)', unique( bases( cb1, : ) ) ) % Test that

% This is expected to find bases belonging to the second circle
cb2 = find_bases2d( edges, verts, bases, d/2, 0, r*1.1 );
assertEquals( (7:12)', unique( bases( cb2, : ) ) ) % Test that
