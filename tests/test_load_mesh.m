function test_load_mesh

fname = '../tests/mesh.txt';
[ verts, edges, bases, epsout, epsin, conductors ] = load_mesh( fname );

assertEquals( [ 2, 1 ], size(conductors) );
assertEquals( transpose( 0:3 ), conductors{1} );
assertEquals( transpose( 4:7 ), conductors{2} );


