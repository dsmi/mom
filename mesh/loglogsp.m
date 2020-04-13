function v = loglogsp( x0, x1, nx )
% v = loglogsp( x0, x1, nx )
%
% Generates a low vector of points with logarithmically increasing steps
% in the first half and logarithmically decreasing in the second

lv = logspace( 1, 2, floor( nx/2 ) + 1 ) - 100;
lw = -logspace( 2, 1, floor( nx/2 ) + 1 );

r = rem( nx, 2 );
vs = [ lv lw(3-r:end) - lw(2-r) ];

v = x0 + ( vs - vs(1) ) * ( x1 - x0 ) / ( vs(end) - vs(1) );
