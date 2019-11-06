function M = mkmommat2l(edges, verts, bases, fintg)
% M = mkmommat2l(edges, verts, bases, fintg)
%
% Computes the moment matrix. The matrix size in NxN, where N is the number
% of the basis functions. The entry M(i,j) corresponds to the integral of
% the Green's function colvolved with base j tested against basis function i.
%
% This one is different from mkmommat2 in that it works with pwl linear basis
% functions which are defined over a pair of edges adjacent to a vertex.
%
%  Params:
%    edges  - num_of_edges-by-2 matrix of the indices of the edge endpoint
%             vertices
%    verts  - num_of_verts-by-2 matrix of the coordinates of the vertices.
%    bases  - num_of_bases-by-2 pair of edges forming a basis function
%    fintg  - handle of the function which evaluates inetgral of the 
%             Green's function. The function should accept two
%             parameters: source and observation segments.
%

N = size( bases, 1 );

M = zeros(N,N);

% To be used as the linear distribution beginning and end values
u0 = zeros(N,1);
u1 = ones(N,1);

% Beginning and end vertices of the edge
ev0 = verts( edges(:,1), : ); 
ev1 = verts( edges(:,2), : );

% Incoming and outgoing edges of a basis function
basein = bases(:,1);
baseout = bases(:,2);

% Matrix is populated by column
for col = 1:N
    
    % Incoming and outgoing source edge, same one for the entire column
    sei0 = repmat( ev0( basein(col),  : ), N, 1 );
    sei1 = repmat( ev1( basein(col),  : ), N, 1 );
    seo0 = repmat( ev0( baseout(col), : ), N, 1 );
    seo1 = repmat( ev1( baseout(col), : ), N, 1 );
    
    % in-to-in, in-to-out, out-to-in, out-to-out
    vii = fintg( sei0, sei1, u0, u1, ev0(basein, :), ev1(basein, :), u0, u1 );
    vio = fintg( sei0, sei1, u0, u1, ev0(baseout,:), ev1(baseout,:), u1, u0 );
    voi = fintg( seo0, seo1, u1, u0, ev0(basein, :), ev1(basein, :), u0, u1 );
    voo = fintg( seo0, seo1, u1, u0, ev0(baseout,:), ev1(baseout,:), u1, u0 );

    M(:,col) = vii + vio + voi + voo;
end
