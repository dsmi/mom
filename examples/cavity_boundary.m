%
% Rectangular cavity 81x61 with two ports on the short sides, boundary
% condition used to improve accuracy
%

addpath( genpath( [ pwd, '/..' ] ) )
addpath( genpath('../../nodal') )

% The geometry
inch2meter = 2.54e-2;
mil2meter = 1.0e-3*inch2meter;
d =  8.0*mil2meter;  % plane-to-plane separation
l = 81.0*mil2meter;  % cavity length
w = 61.0*mil2meter;  % cavity width
fw = 1.0*mil2meter;  % feed lines/port width (one element is w/ny)
a = 1.0*mil2meter;   % metal thickness, used for boundary L and C calculation

% Mesh settings
nx = 80;
ny = 61;

% Dielectric params
lt = 0.02;
er0 = 4.3;
fr = 1e9;

% metal conductivity
sigma = 5.8e7;

% Rectangular cavity, l-by-w
[ e, v ] = mkrect2d(l, w, nx, ny);

% Two ports
c1 = find_edges2d(e, v, -l/2, 0, fw*0.5001);
c2 = find_edges2d(e, v, l/2, 0, fw*0.5001);
%ports = { c1' c2' };
N = size( e, 1 ); % number of boundary segments and ports
ports = mat2cell( 1:N, 1, ones( N, 1 ) ); % ports for the bem solver

%% plotmesh2d( e, v, ports, 1 );

% Frequency-independent circuit solver setup
np = 2;  % number of the ports
nb = N*3;  % number of the branches
branches = [ ones( N, 1 )   (2:N+1)'   ; ...
             (e(:,1)+N+1)   (2:N+1)'   ; ...
             (e(:,2)+N+1)   (2:N+1)' ];

% Branch voltages and currents
W = zeros( nb, np );
K = zeros( nb, np );

% Currents applied to ports
K( c1, 1 ) = 1;
K( c2, 2 ) = 1;

%% % Boundary admittance vs frequency
%% [ fbhz, Sbf ] = SXPParse( 'cavity_edge_8.y1p' );
%% Ybf = s2y( Sbf ) ./ 50; % 50 is the Z0 in cavity_edge_8.y1p
%% fb = fbhz*2*pi;

Y1 = []; % Simulated admittance

% angular frequencies
%freqs = 1e9*2*pi;
freqs = linspace( 1e5, 4e10, 200 )*2*pi;

for freq = freqs,

    % Parameters of the plane.
    er = debye(er0, lt, fr, freq/(2*pi));
    Yplane = j*freq*eps0*er/d;
    Zs = 2*sqrt(j*freq*mu0/sigma); % surface impedance
    Zplane = Zs + j*freq*mu0*d;
    k = sqrt(-Yplane*Zplane);

    fidx = find( freq == freqs );
    nfreqs = length( freqs );
    wavelen = 2*pi/k;
    fprintf( 'Solving for frequency %.8e, %i/%i, wlen %.8e\n', freq, fidx, nfreqs, wavelen );

    fintgsl = @(rsrc,robs) intg_helmsl2d(k,rsrc,robs);
    fintgdl = @(rsrc,robs) intg_helmdl2d(k,rsrc,robs);

    Z = extractz2(e, v, ports, fintgsl, fintgdl);
    Yc = inv(Z*Zplane); % Cavity y-matrix

    % Circuit Y-matrix
    Y = zeros( nb, nb );
    Y( 1:N, 1:N ) = Yc;

    % Find and add the boundary condition

    % Pre-calculated one -- do not use for now
    %% [ mv, midx ] = min( abs( fb - freq ) );
    %% Ybs = Ybf( midx ) * mil2meter; % per-boundary-segment
    
    C = pi*er*eps0/log( d / ( a/2 ) ) % mutual C of two wires
    Ybs = j*freq*C/2*mil2meter; % per-boundary-segment
    Yb = diag( repmat( Ybs, N, 1 ) );
    Yb( c1, c1 ) = 0;
    Yb( c2, c2 ) = 0;
    Y(1:N,1:N) = Y(1:N,1:N) + Yb;

    L = mu0/pi*log( d / ( a/2 ) )
    Zb = j*freq*L*2*mil2meter*0.5; % per-half-boundary-segment
    Y(N+1:end,N+1:end) = Y(N+1:end,N+1:end) + diag( inv(Zb)*ones( N*2, 1 ) );

    Y(N+1:end,N+1:end) = Y(N+1:end,N+1:end) + diag( 1e-10*ones( N*2, 1 ) );

    [ F, V, I ] = solve( branches, Y, W, K );

    % Cavity Y with boundary conditions
    Y = inv( -V( [ c1 c2 ], : ) );

    Y1 = cat(3, Y1, Y);
end

tswrite( 'cavity_mom_81_61_boundary.y2p', freqs/(2*pi), Y1, 'Y', 50 );
