function test_mkbases2d
% Tests mkbases2d
%

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

% Line, dangling vertices at the ends
%
edges = ...
[ 1   2 ;...
  2   3 ;...
  3   4 ];


b0 = [  1   2 ; 2 3 ];

b = mkbases2d( edges );

assertEquals( b0, b );
