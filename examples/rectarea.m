%
% Geometry matches rectarea.hyp, to compare the NaivePdn results against.
%

addpath(genpath([ pwd, '/..' ]));

% The geometry
inch2meter = 2.54e-2;
mil2meter = 1.0e-3*inch2meter;
d = (10.0 + 1.35 + 10.0)*mil2meter; % plane-to-plane separation
l = 1.8*inch2meter;  % cavity length -- x size
w = 0.8*inch2meter;  % cavity width -- y size

% Port radius and locations
x1 = -l/2 + 0.1*inch2meter;
y1 = 0;
x2 = l/2 - 0.1*inch2meter;
y2 = 0;
r = 0.012/2*inch2meter;

% Mesh settings
ny = 40;
nx = round( ny*l/w );
nr = 10;

% Dielectric params
lt = 0.02;
er0 = 4.3;
fr = 1e9;

% metal conductivity
sigma = 5.8e7;

% Rectangular cavity, l-by-w
[ e0, v0 ] = mkrect2d(l, w, nx, ny);

% first port/via
[ e1, v1 ] = mkcir2d(r, nr);
v1 = v1 + repmat( [ x1 y1 ], size( v1, 1 ), 1 ); % move to position
e1(:, [ 1 2 ]) = e1(:, [ 2 1 ]);

% second port/via
[ e2, v2 ] = mkcir2d(r, nr);
v2 = v2 + repmat( [ x2 y2 ], size( v2, 1 ), 1 ); % move to position
e2(:, [ 1 2 ]) = e2(:, [ 2 1 ]);

% Merge all together
e = [ e0; e1 + size(v0, 1); e2 + size(v0, 1) + size(v1, 1) ];
v = [ v0; v1; v2 ];

% Ports
c1 = find_edges2d( e, v, x1, y1, r*1.001 ); % first via
c2 = find_edges2d( e, v, x2, y2, r*1.001 ); % second via
ports = { c1' c2' };

%% plotmesh2d(e,v,ports,1);

Y1 = []; % Simulated admittance

% angular frequencies
%freqs = 1e9*2*pi;
freqs = linspace(1e7, 1e11, 1001)*2*pi;

for freq = freqs,
    freq
    % Parameters of the plane.
    er = debye(er0, lt, fr, freq/(2*pi));
    Yplane = j*freq*eps0*er/d;
    Zs = 2*sqrt(j*freq*mu0/sigma); % surface impedance
    Zplane = Zs + j*freq*mu0*d;
    k = sqrt(-Yplane*Zplane);

    fintgsl = @(rsrc,robs) intg_helmsl2d(k,rsrc,robs);
    fintgdl = @(rsrc,robs) intg_helmdl2d(k,rsrc,robs);

    Z = extractz2(e, v, ports, fintgsl, fintgdl);
    Y = inv(Z*Zplane);
    Y1 = cat(3, Y1, Y);
end

tswrite('rectarea_mom.y2p', freqs/(2*pi), Y1, 'Y', 50);
