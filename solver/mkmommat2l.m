function M = mkmommat2l(bases, fintg)
% M = mkmommat2l(bases, fintg)
%
% Computes the moment matrix. The matrix size in NxN, where N is the number
% of the basis functions. The entry M(i,j) corresponds to the integral of
% the Green's function colvolved with base j tested against basis function i.
%
% This one is different from mkmommat2 in that it works with pwl linear basis
% functions which are defined over a pair of edges adjacent to a vertex.
%
%  Params:
%    bases  - num_of_bases-by-2 pair of edges forming a basis function
%    fintg  - handle of the function which evaluates inetgral of the 
%             Green's function. The function should accept six parameters:
%             source edge, source end weights, observation edge, observation
%             end weights.
%

N = size( bases, 1 );

M = zeros( N, N );

% To be used as the linear distribution beginning and end values
u0 = zeros(N,1);
u1 = ones(N,1);

% Incoming and outgoing edges of a basis function
basein = bases(:,1);
baseout = bases(:,2);

% Matrix is populated by column
for col = 1:N
    
    % Incoming and outgoing source edge, same one for the entire column
    baseincol  = repmat( basein(col) , N, 1 );
    baseoutcol = repmat( baseout(col), N, 1 );
    
    % in-to-in, in-to-out, out-to-in, out-to-out
    vii = fintg( baseincol,  u0, u1, basein,  u0, u1 );
    vio = fintg( baseincol,  u0, u1, baseout, u1, u0 );
    voi = fintg( baseoutcol, u1, u0, basein,  u0, u1 );
    voo = fintg( baseoutcol, u1, u0, baseout, u1, u0 );

    M(:,col) = vii + vio + voi + voo;
end
