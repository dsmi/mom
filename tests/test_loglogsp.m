function test_loglogsp
%

v1 = loglogsp( -1, 1, 2 );
assertEquals( -1, v1(1), 1e-15 );
assertEquals(  1, v1(2), 1e-15 );

v2 = loglogsp( -1, 1, 3 );
assertEquals( -1, v2(1), 1e-15 );
assertEquals(  0, v2(2), 1e-15 );
assertEquals(  1, v2(3), 1e-15 );

v3 = loglogsp( 0, 1, 4 );
assertEquals(  0, v3(1), 1e-15 );
assertEquals(  1.937129433613966e-01, v3(2), 1e-15 );
assertEquals(  8.062870566386034e-01, v3(3), 1e-15 );
assertEquals(  1, v3(4), 1e-15 );

v4 = loglogsp( -1, 1, 5 );
assertEquals(  -1, v4(1), 1e-15 );
assertEquals(  -7.597469266479578e-01, v4(2), 1e-15 );
assertEquals(  0, v4(3), 1e-15 );
assertEquals(  7.597469266479578e-01, v4(4), 1e-15 );
assertEquals(  1, v4(5), 1e-15 );

