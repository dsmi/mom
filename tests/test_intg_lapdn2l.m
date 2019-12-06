function test_intg_lapdn2l
% 
%addpath(genpath([ pwd, '/..' ]));

% The testing geometry - a cicrle
N = 9; % number of segments
[ edges, verts ] = mkcir2d(1, N);
ra = verts(edges(:,1),:);
rb = verts(edges(:,2),:);
ua = reshape( linspace( -3, 4, N ), [], 1 );
ub = reshape( linspace( 2, -1, N ), [], 1 );

[ edges2, verts2 ] = mkcir2d(1.7, N);
verts2 = verts2 * [ 1 0 ; 0 -1 ]; % change sign of y
rc = verts2(edges2(:,1),:);
rd = verts2(edges2(:,2),:);
uc = reshape( linspace( -1, 2, N ), [], 1 );
ud = reshape( linspace( 2, -2, N ), [], 1 );


v = intg_lapdn2l( ra, rb, ua, ub, rc, rd, uc, ud );
vq = intg_lapdn2lq( ra, rb, ua, ub, rc, rd, uc, ud );

assertEquals(vq, v, 1e-6)

% Edges lying on the same plane
ra   = [ 2  1 ];
rb   = [ 3  1 ];
ua   = [   2  ];
ub   = [  -1  ];

rc   = [  -0  1 ];
rd   = [   1  1 ];
uc   = [ -0.5 ];
ud   = [  1.5 ];

v = intg_lapdn2l( ra, rb, ua, ub, rc, rd, uc, ud );
vq = intg_lapdn2lq( ra, rb, ua, ub, rc, rd, uc, ud );

assertEquals(vq, v, 1e-15)


% And the singular edge test
ra   = [ 2  1 ];
rb   = [ 3  1 ];
ua   = [ 1 ];
ub   = [ 1 ];

rc   = [ 2  1 ];
rd   = [ 3  1 ];
uc   = [ 1 ];
ud   = [ 1 ];

v = intg_lapdn2l( ra, rb, ua, ub, rc, rd, uc, ud );

assertEquals(0, v)

