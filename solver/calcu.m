function upt = calcu(edges, verts, intgsl, intgdl, u, q, rpt)
% upt = calcu(edges, verts, intgsl, intgdl, u, q, rpt)
%
%  After the boundary values (both potential and flux) are found, the values
%  at arbitrary point of the domain can be computed. This function computes
%  the potential at a number of points.
%  Inputs:
%   edges  - num_of_edges-by-2 matrix of the indices of the edge endpoint
%             vertices
%   verts  - num_of_verts-by-2 matrix of the coordinates of the vertices.
%   ports  - cell array of vectors, edges of the ports.
%   intgsl - handle of a function which evaluates integrals of single layer
%            potential, is used to call mkmommat2.
%   intgdl - handle of a function which evaluates integrals of double layer
%            potential, is used to call mkmommat2.
%   u, q   - num_of_ports-by-num_of_edges boundary potentials and fluxes.
%            As computed by (for example) extractz2
%   rpt    - n-by-2, the evaluation points
% Outputs:
%   upt    - The resulting potentials at rpt
%

% Number of points
npt = size(rpt, 1);

% Number of edges in the geo
N = size(edges,1);

m = repmat((1:npt)', 1, N);
n = repmat((1:N), npt, 1);

m = m(:);
n = n(:);
rsrc = cat(3, verts(edges(n,1),:), verts(edges(n,2),:));
robs = cat(3, rpt(m,:), rpt(m,:)); % hack - we know that the integration
                                   % routines just find the center


G = zeros(npt,N);
G(:) = intgsl(rsrc, robs);

H = zeros(npt,N);
H(:) = intgdl(rsrc, robs);

upt = G*q-H*u;
