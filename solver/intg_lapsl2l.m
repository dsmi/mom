function v = intg_lapsl2l( ra, rb, ua, ub, rc, rd, uc, ud )
% v = intg_lapsl2l( ra, rb, ua, ub, rc, rd, uc, ud )
%
% Evaluates integral of the potential due to the charge linearly distributed
% over the sorce edge along the obvervation edge, weighted by the linear
% funtion:
%  v = \int_obs f(u) \int_src f(v) -1/(2*pi)*ln(|v-u|/c) dv du
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

% This pre-allocates the output array
v = 0 * ua;

% Quadrature is used to integrate over the observation segments
[qx,qw] = GLTable(7);
%% [qx,qw] = GLNodeWt(21);

for xw = [ reshape( qx, 1, [] ) ; reshape( qw, 1, [] ) ]

    robs = rc + ( rd - rc ) * ( xw(1) * 0.5 + 0.5 );
    uobs = uc + ( ud - uc ) * ( xw(1) * 0.5 + 0.5 );

    roa = ra - robs;
    rob = rb - robs;

    % We use a coordinate system with zero at the observation point and x
    % axis parallel to the source edge. Here we translate to this system.
    h = sum( roa .* n, 2 );
    a = sum( roa .* t, 2 );
    b = sum( rob .* t, 2 );

    la2 = h.^2 + a.^2;
    lb2 = h.^2 + b.^2;

    c = ( b.*ua - a.*ub ) ./ l;
    g = ( ub - ua ) ./ l;

    % This is to avoid division by zero. h can be negative, but it is unlikely
    % to ever be equal to -1e-300
    hnz = h + 1e-300;

    vq = 1/4*( log(lb2).*( 2*c.*b + g.*lb2 ) - log(la2).*( 2*c.*a + g.*la2 ) ...
          - 4*l.*c - g.*( b.*b - a.*a ) + 4*c.*h.*( atan(b./hnz) - atan(a./hnz) ) );

    v = v + xw(2) .* uobs .* obsl./2 .* vq;

end

% Intergration over X - Constant
%v = l .* ( log( l./2 ) - 1 );
% fintg = @( x, y ) 1/2 * ( x - y ) * ( log( ( x - y )^2 ) - 2 );
%v = fintg( l, y ) - fintg( 0, y );
%% v = 1/2 * ( ( l - y ) * ( log( ( l - y )^2 ) - 2 ) + y * ( log( ( -y )^2 ) - 2 ) );

% Intergration over X - Linear
%v = -1/2 * l .* ( -log( l ) + 1 + log( 2 ) );
%% fintg = @( x, y ) 1/4 * (x-y) .* ( (x+y) .* log((x-y).^2) - x - 3*y )
%% v = ( fintg( l, y ) - fintg( 0, y ) ) ./ l;
%% v = 1./(4*l).* ( (l-y) .* ( (l+y) .* log( (l-y).^2) - l - 3*y ) - (-y) .* ( y.*log((-y).^2) - 3*y ) );

% Intergration over X - Linear with the given end values:
% u(x) = u1 + (u2-u1)*x/l
%% k = (ub-ua) ./ (4*l);
%% v =  ( ua./2.*(l-y) + k.*(l-y).*(l+y) ).*log( (l-y).^2 ) ...
%%     + ( ua./2.*y  + k.*y.*y ) .* log( ( -y ).^2 ) ...
%%     - 2*k.*l.*y ...
%%     - (3*ua+ub).*l./4 ;

%% % Integration over y, constant:
%% k = (ub-ua) ./ (4*l);
%% v = 1/36 * l.*l .* ( 3*log( l.*l ) .* ( 3*ua + 8*k.*l ) - 9*ua - 28*k.*l ) ...
%%     + 1/36 * l.*l .* ( 3*log( l.*l ).*( 4*k.*l+3*ua ) - 8*k.*l - 9*ua ) ...
%%     - k.*l.^3 ...
%%     - (3*ua+ub).*l.*l./4 ;

%% % Linear
%% k = (ub-ua) ./ (4*l);
%% v = 1/72 * l.^3 .* ( 6*log( l.*l ) .* ( ua + 3*k.*l ) - 10*ua - 33*k.*l ) ...
%%     + 1/72 * l.^3 .* ( 6*log( l.*l ) .* ( 3*k.*l + 2*ua ) - 9 * k.*l - 8*ua ) ...
%%     - 2*k.*l.^4/3 ...
%%     - (3*ua+ub).*l.*l.*l./8 ;
%% v = v./l;

% Over Y linear with the given end values:
% u(x) = u1 + (u2-u1)*x/l

% Singular terms here
vs = l.^2 .* ( 1/8*log( l.^2 ).*(ud+uc).*(ua+ub) ...
        - (7/16*ua+5/16*ub).*uc - (7/16*ub+5/16*ua).*ud );

% Find singular pairs
sidx = find( all( abs(ra-rc)<1e-12 & abs( abs(rb-rd)<1e-12 ), 2) );

v(sidx) = vs(sidx);

%% % Zero zero-len edges
%% v( find( all( ra == rb, 2 ) | all( rc == rd, 2 ) ) ) = 0;
