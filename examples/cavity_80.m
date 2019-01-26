%
% Rectangular cavity with two ports on the narrow side, to be combined
% with the stripline models -- power aware si experiments.
%

addpath(genpath([ pwd, '/..' ]));

% The geometry
inch2meter = 2.54e-2;
mil2meter = 1.0e-3*inch2meter;
d = (5.0)*mil2meter; % plane-to-plane separation
w = 80.0*mil2meter;  % cavity width
l = 80.0*mil2meter;  % cavity length

% Mesh settings
nx = 80;
ny = 80;
pw = 2.0*mil2meter; % port width; (one element is w/ny)

% Dielectric params
lt = 0.02;
er0 = 4.3;
fr = 1e9;

% metal conductivity
sigma = 5.8e7;

% Rectangular cavity, l-by-w
[ e, v ] = mkrect2d(l, w, nx, ny);

% Two ports in the middle of the opposite sides
c1 = find_edges2d(e, v, -l/2, 0, pw*0.5001);
c2 = find_edges2d(e, v, l/2, 0, pw*0.5001);
ports = { c1' c2' };

%% plotmesh2d(e,v,ports,1);

Y1 = []; % Simulated admittance

% angular frequencies
%freqs = 1e9*2*pi;
freqs = linspace(1e7, 4e10, 500)*2*pi;

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

tswrite('cavity_80_port2.y2p', freqs/(2*pi), Y1, 'Y', 50);
