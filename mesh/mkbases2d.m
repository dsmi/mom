function b = mkbases2d(e)
% Defines 'bases' or support of the basis functions. A base is a vertex
% of the mesh with one incoming and one outgoing edges.
% Output is N-by-2 array, first column is index of the incoming edge (second
% vertex of this edge is shared) and second column is the outgoing vertex
% index.
% The current implementation assumes that each vertex has a base
% associated with it.

% Incoming and outgoing edges of a vertex (assuming no vertices
% with more than two edges connected)
ei = eo = zeros( max( e(:) ), 1 );
ei(e(:,2)) = 1:size(e,1);
eo(e(:,1)) = 1:size(e,1);

% Indices of the vertices with two edges connected.
v = find( ei > 0 & eo > 0);

% Finally, bases!
b = [ v*0 v*0 ];
b(:,1) = ei(v);
b(:,2) = eo(v);
