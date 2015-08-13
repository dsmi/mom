function img = mkimage(edges, verts, intgsl, intgdl, u, q, nx, ny)
% img = mkimage(edges, verts, intgsl, intgdl, u, q, nx, ny)
%
% Builds an image of the potential - calculates the potential at 2d square
% array of points. Calls calcu to calculate the potential values.
% Inputs:
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
%   nx, ny - dimensions of the image
% Outputs:
%   upt    - The resulting image, ny-by-nx array.

% Obtain the geometry bounding box
boxmin = min(verts);
boxmax = max(verts);
boxsize = boxmax - boxmin;

% Scale to translate from pixels to geo
scale = boxsize./[ nx-1 ny-1 ];

% Pixel indices, one-based
[ px, py ] = meshgrid(1:nx, 1:ny);

% Populate the image by columns to avoid memory overflow
img = zeros(ny,nx);

for l=1:nx

    % Pixel coordinates
    rpix = [ px(:,l)-1 py(:,l)-1 ].*repmat(scale, ny, 1) + repmat(boxmin, ny, 1);

    img(:,l) = calcu(edges, verts, intgsl, intgdl, u, q, rpix);

end
