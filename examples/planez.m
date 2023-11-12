%
% Surface impedance of a plane of the given thickness.
%

addpath(genpath([ pwd, '/..' ]));

% rod params
ccopper = 5.8e7;
sigma = ccopper; % conductivity -- copper
w = 3e-3; % width
t = 1e-4; % thickness
nw = 200;
nt = 3; 

[ edges, verts ] = mkrect2d( w, t, nw, nt );

c1 = find_eol2d( edges, verts, 0, 1, t/2);
c2 = find_eol2d( edges, verts, 0, -1, t/2);

ports = { c1' };%c2' };

%% plotmesh2d( edges, verts, ports, 1 );

freqs = 2*pi*1e6;%linspace(5e8,5e9,3);

for freq = freqs,

    k = (1-j)*sqrt( freq*mu0*sigma/2 ) % wavenumber
    k = freq * sqrt(mu0 * (eps0 - j*sigma/freq))
    
    nu = sqrt( freq*mu0/sigma )*exp( j*pi/4 ) % impedance
    nu = sqrt( mu0/(eps0 - j*sigma/freq) ) % impedance
    delta = sqrt(2./(freq*mu0*sigma)) % skin depth

    fintgsl = @(rsrc,robs) intg_helmsl2d(k,rsrc,robs);
    fintgdl = @(rsrc,robs) intg_helmdl2d(k,rsrc,robs);

    % Condutor medium
    z = j*freq*mu0;
    y = j*freq*eps0 + ccopper;

    Z = extractz2( edges, verts, ports, fintgsl, fintgdl )*w*z

    k = sqrt( -z.*y )
    nu = sqrt( z./y )

    Z11 = nu.*(1+exp(-j*k*2*t))./(1-exp(-j*k*2*t))
    Z21 = nu*2./(exp(j*k*t)-exp(-j*k*t))
    
    
end

