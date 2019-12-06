function v = intg_lapdn2lq( ra, rb, ua, ub, rc, rd, uc, ud )
% v = intg_lapdl2lq( ra, rb, ua, ub, rc, rd, uc, ud )
%
% Double layer field calculated with quadrature, to be used
% for testing.
% 
%  ra, rb - endpoints of the source segments, N-by-2
%  ua, ub - f(v) at the ends of the source segments, N-by-1
%  rc, rd - endpoints of the observation segments, N-by-2
%  uc, ud - f(u) at the ends of the observation segments, N-by-1
% 
%  v - the resulting integrals, N-by-1
%

% Number of edges
N = size(ra,1);

% Source edge vectors. N-by-2 array.
se = rb - ra;

% Source edge lengths. Column vector of length N.
sl = sqrt( sum( se.^2, 2 ) );

% Source edge tangentials - normalized edge vectors.
st = se ./ sl(:,ones(1,2));

% Observation edge vectors. N-by-2 array.
oe = rd - rc;

% Observation edge lengths. Column vector of length N.
ol = sqrt( sum( oe.^2, 2 ) );

% Observation edge tangentials - normalized edge vectors.
ot = oe ./ ol(:,ones(1,2));

% Observation outer normal. The outer polygon is CCW, the holes are CW,
% thus to get the outer normals we need to rotate the edge tangentials
% 90 degrees clockwise.
n = [ ot(:,2) -ot(:,1) ];

% Number of quadrature points used.
qN=11;

[qX,qW] = GLNodeWt(qN);

% This pre-allocates the output array
v = 0 * ua;

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

            % Source-to-observation vector
            r = ro - rs;

            % Squared source-to-observation distance
            R2 = sum( r.^2, 2 );

            % n*r/R^2 
            k = us * uo * sxw(2)*sl(ei)/2 * oxw(2)*ol(ei)/2;
            v( ei ) = v( ei ) + k * sum( n(ei,:).*r, 2 ) ./ R2;
            
        end
    end
end
