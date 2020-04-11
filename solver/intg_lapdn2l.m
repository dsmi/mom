function v = intg_lapdn2l( ra, rb, ua, ub, rc, rd, uc, ud )
% v = intg_lapdn2l( ra, rb, ua, ub, rc, rd, uc, ud )
%
% Evaluates the following double integral:
%  v = \int_obs f(u) \int_src f(v) -1/(2*pi) (r*n) / R^2  dv du
%
% Where f(u) and f(v) are the weighting functions, r is the
% source-to-observation vector, R = |r|.
% Depending on the need/context, it ca be viewed as:
%  (a) Electric field due to the charge on the source edge, dotted with
%      the observation normal (can be used for enforcing continuity of
%      the electric field at the dielectric-dielectric interface)
%  (b) Potential at the source due to the double electric layer at
%      observation.
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

% Edge vectors. N-by-2 array.
edges = rb - ra;

% Edge lengths. Column vector of length N.
l = sqrt( sum(edges.^2,2) );

% Edge tangentials - normalized edge vectors.
t = edges ./ l(:,ones(1,2));

% Edge normals -- tangentials rotated counterclockwise.
n = [ -t(:,2) t(:,1) ];

% Observation edge lengths. Column vector of length N.
obsl = sqrt( sum( ( rd - rc ).^2, 2 ) );

% Observation edge vectors. N-by-2 array.
oe = rd - rc;

% Observation edge lengths. Column vector of length N.
ol = sqrt( sum( oe.^2, 2 ) );

% Observation edge tangentials - normalized edge vectors.
ot = oe ./ ol(:,ones(1,2));

% Observation inner normal. Notice that if the outer polygon is CCW then
% this is the inner normal (sign of the result might need to be changed)
on = [ -ot(:,2) ot(:,1) ];

% This pre-allocates the output array
v = 0 * ua;

% Quadrature is used to integrate over the observation segments
[qx,qw] = GLTable(7);
%% [qx,qw] = GLNodeWt(31);

for xw = [ reshape( qx, 1, [] ) ; reshape( qw, 1, [] ) ]

    robs = rc + ( rd - rc ) * ( xw(1) * 0.5 + 0.5 );
    uobs = uc + ( ud - uc ) * ( xw(1) * 0.5 + 0.5 );

    % observation-to-a and observation-to-b vectors
    roa = ra - robs;
    rob = rb - robs;

    % Squared distances from a and b to observation
    R2a = sum( roa.^2, 2 );
    R2b = sum( rob.^2, 2 );

    % We use a coordinate system with zero at the observation point and x
    % axis parallel to the source edge. Here we translate to this system.
    h = sum( roa .* n, 2 );
    a = sum( roa .* t, 2 );
    b = sum( rob .* t, 2 );

    g = ( ub - ua ) ./ l;

    % ingegral of x / R^2
    fx = 1/2 * ( log( R2b ) - log( R2a ) );

    % This is to avoid division by zero. h can be negative, but it is unlikely
    % to ever be equal to -1e-300
    hnz = h + 1e-300;

    % ingegral of y / R^2
    fy = ( atan(b./hnz) - atan(a./hnz) );
    
    vx = ( ua - g.*a ) .* fx + g.*( ( b - a ) - h .* fy );
    vy = ( ua - g.*a ) .* fy + g.*h.*fx;
    
    % vx is along t and vy is along n, project both on
    % the observation normal (on)
    vq = vx .* dot( on, t, 2 ) + vy .* dot( on, n, 2 );

    v = v + xw(2) .* uobs .* obsl./2 .* vq;

end

% Find singular pairs
sidx = find( all( abs(ra-rc)<1e-12 & abs( abs(rb-rd)<1e-12 ), 2) );

% Integral of the normal derivative is zero for the self terms -- due
% to the symmetry.
v(sidx) = 0*v(sidx);

%% % Zero zero-len edges
%% v( find( all( ra == rb, 2 ) | all( rc == rd, 2 ) ) ) = 0;

