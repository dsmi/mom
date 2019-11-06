%function test_intg_lapsl2l
%
% Test intg_lapsl2l
%

ra   = [ -0.5  4  ; -0.5  4 ; -0.5  4  ; -0.5  4 ; 1 2 ];
rb   = [  0.5  6  ;  0.5  6 ;  0.5  6  ;  0.5  6 ; 2 3 ];
ua   = [  1  ; 0 ; 2 ; -1 ; 0 ];
ub   = [  0  ; 1 ; 1 ; -2 ; 1 ];

rc   = [  -0.5  -2 ;  -0.5  -2 ; -0.5  -2 ;  -0.5  -2 ; 2 3 ];
rd   = [   0.5  -3 ;  0.5  -3  ;  0.5  -3 ;  0.5  -3  ; 3 4 ];
uc   = [ 1 ; 0 ; -2 ; 2 ; 1 ];
ud   = [ 0 ; 1 ;  1 ; 1 ; 0 ];

v = intg_lapsl2l( ra, rb, ua, ub, rc, rd, uc, ud );

% This pre-allocates the test values array
vt = 0 * ua;

% Number of quadrature points used.
qN=10;

[qX,qW] = GLNodeWt(qN);

% Number of edges
N = size(ra,1);

% Source edge lengths. Column vector of length N.
sl = sqrt( sum( ( rb - ra ).^2, 2 ) );

% Observation edge lengths. Column vector of length N.
ol = sqrt( sum( ( rd - rc ).^2, 2 ) );

% Edges loop
for ei = 1:N

    % Source edge quadrature
    for sxw = [ reshape( qX, 1, [] ) ; reshape( qW, 1, [] ) ]
        
        % Observation edge quadrature
        for oxw = [ reshape( qX, 1, [] ) ; reshape( qW, 1, [] ) ]

            % Source point and linear weight value
            rs = ra(ei,:) + ( rb(ei,:) - ra(ei,:) ) * ( sxw(1) * 0.5 + 0.5 );
            us = ua(ei)   + ( ub(ei)   - ua(ei,:) ) * ( sxw(1) * 0.5 + 0.5 );

            % Observation point and linear weight value
            ro = rc(ei,:) + ( rd(ei,:) - rc(ei,:) ) * ( oxw(1) * 0.5 + 0.5 );
            uo = uc(ei)   + ( ud(ei)   - uc(ei,:) ) * ( oxw(1) * 0.5 + 0.5 );

            % Source-to-observation distance
            R = sqrt( sum( ( ro - rs ).^2, 2 ) );

            k = us * uo * sxw(2)*sl(ei)/2 * oxw(2)*ol(ei)/2;
            vt( ei ) = vt( ei ) + k * log( R );
            
        end
    end
end

% Fifth pair are close by, so the accuracy is lower
assertEquals(vt(1:4), v(1:4), 1e-13);
assertEquals(vt(5), v(5), 2e-4);

% The non-self-term calculation procedure, which involves the quadrature,
% works for the self-terms as well. So in order to test the self-term
% calculation we first call the integral evaluation
% function for the exact same source and observation edges -- so the self-term
% calculation is used, and then we offset the observations a little bit so the
% non-self-term is used, and compare the results.

ra   = [ -0.5  4 ; -0.5  4 ; -0.5  4 ];
rb   = [  0.5  6 ;  0.5  6 ;  0.5  6 ];
ua   = [  1 ; 0 ; 2 ];
ub   = [  0 ; 1 ; 4 ];

rc   = [ -0.5  4 ; -0.5  4 ; -0.5  4 ];
rd   = [  0.5  6 ;  0.5  6 ;  0.5  6 ];
uc   = [ 0 ; 1 ; 3 ];
ud   = [ 1 ; 0 ; -2 ];

% Self-terms are calculated here
vs = intg_lapsl2l( ra, rb, ua, ub, rc, rd, uc, ud );

% Offset applied here to use the 'regular' calculation with quadrature
vst = intg_lapsl2l( ra, rb, ua, ub, rc + 1e-11, rd + 1e-11, uc, ud );

assertEquals(vst, vs, 1e-3);

% To make sure the different calculation routines are used
assertTrue( ~nnz( vst == vs ) );

