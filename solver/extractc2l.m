function C = extractc2l( edges, verts, bases, epsout, epsin, conductors )
% C = extractc2l( edges, verts, bases, epsout, epsin, conductors )
% 
% Calculates N-by-N capacitance matrix of N-conductor system similar to
% what extractc2 does, but using piecewise-linear basis functions.
%
% Inpnputs:
%   edges      - num_of_edges-by-2 matrix of the indices of the edge endpoint
%                vertices
%   verts      - num_of_verts-by-2 matrix of the coordinates of the vertices.
%   bases      - num_of_bases-by-2 pair of edges forming a basis function
%   epsout     - column vector of lenght num_of_edges, permittivity outside
%                the given edge: the side where the normal of the edge points
%                is considered outside.
%   epsin      - column vector of lenght num_of_edges, permittivity inside.
%                Only makes sense for dielectric-to-dielectric boundary edges.
%   conductors - cell array of vectors, edges of conductors.
% Outputs:
%   C          - the resulting capacitance matrix.
%

% Number of the edges
ne = size( edges, 1 );

% Number of the bases
nb = size( bases, 1 );

% Bases forming the conductor boundaries
conductor_bases = cellfun( @( c ) find( ismember( bases(:,1), c ) ), ...
                           conductors, 'UniformOutput', false );

% Vector of length nb, 1 if it is the conductor base, 0 otherwise
iscnd = zeros( nb,1);
iscnd( cell2mat( conductor_bases ) ) = 1;

% Edge lengths are needed to compute the potential integrals
% and the total charge per each conductor. N-by-1
r1 = verts(edges(:,1),:);
r2 = verts(edges(:,2),:);
edgelen = sqrt( sum( ( r2 - r1 ).^2, 2) );

% Length of the edge pair forming the base
baselen = edgelen( bases(:,1) ) + edgelen( bases(:,2) );

% The matrix equation enforces potentials on the conductor edges and normal derivative
% of the electric field on the dielectric edges. It is:
%    [ P ; E ]*q=[ p ; 0 ]
% where q is the per-length segment charges, and p is potetnials.
% P are the potential coefficients and E are the electric field
% coefficients.

% Beginning and end vertices of the edge
ev0 = verts( edges(:,1), : ); 
ev1 = verts( edges(:,2), : );

% Multiplier to be applied to the integrals
m = (1/eps0)*(-1/(2*pi));

% Integration funtion wrapper so the mkmommat only needs to pass the edges
intgp = @( srcedge, ua, ub, obsedge, uc, ud ) ...
         m * intg_lapsl2l( ev0( srcedge, : ), ev1( srcedge, : ), ua, ub, ...
                           ev0( obsedge, : ), ev1( obsedge, : ), uc, ud );

% P matrix, charge-to-potential
P = mkmommat2l( bases, intgp );

% Wrapper for intg_lapdn2l, sign is changed because intg_lapdn2l assumes
% different edges orientation
intge = @( srcedge, ua, ub, obsedge, uc, ud ) ...
        - m * intg_lapdn2l( ev0( srcedge, : ), ev1( srcedge, : ), ua, ub, ...
                            ev0( obsedge, : ), ev1( obsedge, : ), uc, ud );

% Self-terms of E matrix calculation
intgs = @( srcedge, ua, ub, obsedge, uc, ud ) ...
         intg_lapdn2s( edgelen, epsout, epsin, ...
                       srcedge, ua, ub, obsedge, uc, ud );

% E matrix, charge to normal D (displacement) discontinuity
E = mkmommat2l( bases, intge ) + mkmommat2l( bases, intgs );

% lhs matrix, elements of P/E are used for conductor/dielecteic boundary
A = diag(iscnd)*P + diag(1-iscnd)*E;

% Number of the conductors.
nc = length(conductors);

% Build p matrix - the number of columns corresponds to the number
% of conductors, each column has p=1 for one of the conductors and
% p=0 for all the others.
p = zeros( nb, nc );
[ baseidx bportidx ] = ports2subs( conductor_bases );
p( sub2ind(size(p), baseidx, bportidx) ) = 1/2*baselen( baseidx );

% Find the base charges
qb = A\p;

% A matrix to find segment/edge charges from the base charges
T = 0.5 * ( sparse( bases(:,1), 1:nb, ones( nb, 1 ), ne, nb ) ...
            + sparse( bases(:,2), 1:nb, ones( nb, 1 ), ne, nb ) );

% Edge/segment charges
q = T*qb;

% The chagres found on the conductor boundaries are the total ones,
% convert the total charge to free charge
qf = diag(epsout./eps0)*q;

% Q matrix multiplies the edge charge densities by the edge length
% to get the total charge associated with the given edge, and sums
% the edge charges to get the conductor charges.
% Q(m,n) = length(n) if n belongs to port(m) and = 0 otherwise.
Q = zeros( nc, ne );
[ edgeidx eportidx ] = ports2subs( conductors );
Q( sub2ind(size(Q), eportidx, edgeidx) ) = edgelen( edgeidx );

C = Q*qf;
