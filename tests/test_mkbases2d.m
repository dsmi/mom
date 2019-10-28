function test_mkbases2d
% Tests mkbases2d
%

% Current mkbases2d does not support dangling vertices.
% Give it a simple test though.

edges = ...
[ 3   2 ;...
  2   1 ;...
  1   3 ];


b0 = ...
[  2   3 ;...
   1   2 ;...
   3   1 ];

b = mkbases2d( edges );

assertEquals( b0, b );
