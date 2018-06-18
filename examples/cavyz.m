%
% Y and Z elements of a section of the cavity equivalent circuit.
%

addpath(genpath([ pwd, '/..' ]));

% The geometry
inch2meter = 2.54e-2;
mil2meter = 1.0e-3*inch2meter;
d = (5.0 + 0.675 + 5.0)*mil2meter; % plane-to-plane separation
w = 100.0*mil2meter;  % cavity width
l = 200.0*mil2meter;  % cavity length

% Mesh settings
ny = 50;
nx = ny*2;

% Cell size
a = l/nx
b = w/ny

% Dielectric params
lt = 0.02;
er0 = 4.3;
fr = 1e9;

% metal conductivity
sigma = 5.8e7;

Y1 = [];
Z1 = [];

fhz0 = 0;
fhz1 = 5e10;
%freqs = logspace(log10(fhz0*2*pi),log10(fhz1*2*pi),2000);
freqs = linspace(fhz0*2*pi,fhz1*2*pi,2000);

for freq = freqs,
    % Parameters of the plane.
    er = debye(er0, lt, fr, freq/(2*pi));
    Yplane = j*freq*eps0*er/d;
    Zs = 2*sqrt(j*freq*mu0/sigma); % surface impedance
    Zplane = Zs + j*freq*mu0*d;

    Y1 = cat(3, Y1, Yplane*a*b);
    Z1 = cat(3, Z1, Zplane*a/b*0.5); % each section has four half-Z
end

tswrite('cavy.y1p', freqs/(2*pi), Y1, 'Y', 50)
tswrite('cavz.z1p', freqs/(2*pi), Z1, 'Z', 50)
