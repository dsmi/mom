%
% Rectangular cavity 81x61 with two ports on the short sides.
%

addpath(genpath([ pwd, '/..' ]));

% The geometry
inch2meter = 2.54e-2;
mil2meter = 1.0e-3*inch2meter;
d =  8.0*mil2meter;  % plane-to-plane separation
l = 81.0*mil2meter;  % cavity length
w = 61.0*mil2meter;  % cavity width
fw = 1.0*mil2meter;  % feed lines/port width (one element is w/ny)

% Mesh settings
nx = 81;
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
ports = { c1' c2' };

%% plotmesh2d( e, v, ports, 1 );

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
    Y = inv(Z*Zplane);
    Y1 = cat(3, Y1, Y);
end

tswrite('cavity_mom_81_61.y2p', freqs/(2*pi), Y1, 'Y', 50);
