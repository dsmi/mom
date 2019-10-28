function b = mkbases2d(e)
% Defines 'bases' or support of the basis functions. A base is a vertex
% of the mesh with one incoming and one outgoing edges.
% Output is N-by-2 array, first column is index of the incoming edge (second
% vertex of this edge is shared) and second column is the outgoing vertex
% index.
% The current implementation assumes that each vertex has a base
% associated with it.

% Number of vertices
nv = max(e(:));
    
b1 = b2 = -ones( nv, 1 );

b1( e(:,2) ) = 1:nv;
b2( e(:,1) ) = 1:nv;

b = [ b1, b2 ];
