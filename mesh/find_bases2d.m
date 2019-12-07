function b = find_bases2d( bases, eidx )
% b = find_bases2d( bases, eidx )
%
% Finds the basis funtions (defined over a pair of edges adjecent to a vertex)
% formed by the given set of edges.
% a given point.
%
% Inputs:
%    bases  - num_of_bases-by-2 pair of edges forming a basis function
%    eidx   - indides of the targed edges.
% Outputs:
%    b      - vector of the found base indices

% Sort for ismember
e = sort( eidx );

% And find bases with each edges in the circle
b = find( ismember( bases(:,1), e ) & ismember( bases(:,2), e ) );

    
