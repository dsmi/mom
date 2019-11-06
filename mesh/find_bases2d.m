function b = find_bases2d( edges, verts, bases, x, y, radius )
% b = find_bases2d( edges, verts, bases, x, y, radius )
%
% Finds the basis funtions (defined over a pair of edges adjecent to a vertex)
% with each of the support edges located within a particular radius around
% a given point.
%
% Inputs:
%    edges  - num_of_edges-by-2 matrix of the indices of the edge endpoint
%             vertices
%    verts  - num_of_verts-by-2 matrix of the coordinates of the vertices.
%    bases  - num_of_bases-by-2 pair of edges forming a basis function
%    x, y   - center point of the lookup area
%    radius - radius of the lookup area
% Outputs:
%    b      - vector of the found base indices

% Find conforming edges first, and sort for ismember
e = sort( find_edges2d( edges, verts, x, y, radius ) );

% And find bases with each edges in the circle
b = find( ismember( bases(:,1), e ) & ismember( bases(:,2), e ) );

    
