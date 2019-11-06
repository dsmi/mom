function C = extractc2l(edges, verts, bases, epsd, conductors)
% C = extractc2l(edges, verts, bases, epsd, conductors)
% 
% Calculates N-by-N capacitance matrix of N-conductor system, embedded in
% uniform dielectric, using piecewise-linear basis functions.
%
% Inpnputs:
%   edges      - num_of_edges-by-2 matrix of the indices of the edge endpoint
%                vertices
%   verts      - num_of_verts-by-2 matrix of the coordinates of the vertices.
%   bases      - num_of_bases-by-2 matrix, each row is a pair of the edge
%                indices forming the corresponding basis function.
%                See mkbases2d for additional details.
%   epsd       - permittivity of the dielectric
%   conductors - cell array of vectors, edges of conductors.
% Outputs:
%   C          - the resulting capacitance matrix.
%

% Number of the edges
N = size(edges, 1);

% P are the potential coefficients
P = (1/epsd)*(-1/(2*pi))*mkmommat2l(edges, verts, bases, @intg_lapsl2l);

% lhs matrix = P
A = P;

% Number of the conductors.
nc = length(conductors);

% Edge lengths are needed to compute the potential integrals
% and the total charge per each conductor. N-by-1
r1 = verts(edges(:,1),:);
r2 = verts(edges(:,2),:);
edgelen = sqrt( sum( ( r2 - r1 ).^2, 2) );

% Length of the edge pair forming the base
baselen = edgelen( bases(:,1) ) + edgelen( bases(:,2) );

% Build p matrix - the number of columns corresponds to the number
% of conductors, each column has p=1 for one of the conductors and
% p=0 for all the others.
p = zeros(N,nc);
[ faceidx fportidx ] = ports2subs( conductors );
p( sub2ind(size(p), faceidx, fportidx) ) = 1/2*baselen( faceidx );

% Find the charges
q = A\p;

% The chagres found on the conductor boundaries are the total ones,
% convert the total charge to free charge
qf = q;

% Q matrix multiplies the basis weights by the base support length
% divided by two -- to get the total charge associated with the
% given base, and sums the base charges to get the conductor charges.
% Q(m,n) = length(n)/2 if n belongs to port(m) and = 0 otherwise.
Q = zeros(nc,N);
Q( sub2ind(size(Q), fportidx, faceidx) ) = 1/2*baselen( faceidx );

C = Q*qf;
