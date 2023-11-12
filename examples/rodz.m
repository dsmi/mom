%
% Surface impedance of a round rod.
%

addpath(genpath([ pwd, '/..' ]));

% rod params
sigma = 5.8e7; % conductivity -- copper
a = 2e-5;  % radius
n = 500; % number of edges

[ edges, verts ] = mkcir2d(a, n);

ports = { 1:size( edges, 1 ) };

%% plotmesh2d( edges, verts, ports, 1 );

Z1 = []; % Simulated impedance
Z2 = []; % Analytical

freqs = 2*pi*linspace(1e3,1e9,50);

for freq = freqs,

    freq
    k = (1-j)*sqrt( freq*mu0*sigma/2 ) % wavenumber
    
    nu = sqrt( freq*mu0/sigma )*exp( j*pi/4 ) % impedance
    nu = sqrt( mu0/(eps0 - j*sigma/freq) ) % impedance
    delta = sqrt(2./(freq*mu0*sigma)) % skin depth

    % Condutor medium
    z = j*freq*mu0;
    y = j*freq*eps0 + sigma;

    fintgsl = @(rsrc,robs) intg_helmsl2d(k,rsrc,robs);
    fintgdl = @(rsrc,robs) intg_helmdl2d(k,rsrc,robs);

    Z = extractz2( edges, verts, ports, fintgsl, fintgdl )*z

    R0 = 1/(pi*a*a*sigma);
    kx = ( 1 - j ) ./ delta;
    Zr = R0 .* kx .* a/2 .* besselj( 0, kx.*a ) ./ besselj( 1, kx.*a )
    
    Z1 = [ Z1 Z ];
    Z2 = [ Z2 Zr ];
    
end

plot( freqs, real(Z1), '-r', freqs, real(Z2), '*r', freqs, imag(Z1), '-b', freqs, imag(Z2), '*b' );
