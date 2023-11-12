%
% Surface impedance of a tube
%

addpath(genpath([ pwd, '/..' ]));

% rod params
sigma = 5.8e7; % conductivity -- copper
r1 = 5e-5;  % inner radius
r2 = 6e-5;  % outer radius
n = 500; % number of edges

[ e1, v1 ] = mkcir2d(r1, n);
e1(:, [ 1 2 ]) = e1(:, [ 2 1 ]);

[ e2, v2 ] = mkcir2d(r2, n);

edges = [ e1; e2 + size(v1, 1) ];
verts = [ v1; v2 ];

% join meshes
edges = [ e1 ; e2 + size(v1, 1) ];
verts = [ v1 ; v2 ];

% Inner circle is the port
ports = { 1:size( e1, 1 ) };

%% plotmesh2d( edges, verts, ports, 1 );

Z1 = []; % Simulated impedance
Z2 = []; % Analytical

freqs = 2*pi*linspace(1e5,7e8,5);

for freq = freqs,

    freq

    delta = sqrt(2./(freq*mu0*sigma)); % skin depth

    % Condutor medium
    z = j*freq*mu0;
    y = j*freq*eps0+sigma;
    k = sqrt( -y*z );

    fintgsl = @(rsrc,robs) intg_helmsl2d(k,rsrc,robs);
    fintgdl = @(rsrc,robs) intg_helmdl2d(k,rsrc,robs);

    Z = z*extractz2(edges, verts, ports, fintgsl, fintgdl);

    Z1 = [ Z1 Z ];

    % Analytical solution
    if freq/(2*pi) < 2e8,
        
        M = [ besselh(0,2,k*r1) besselh(0,1,k*r1) ; ...
              besselh(0,2,k*r2) besselh(0,1,k*r2) ];

        dM = [ -besselh(1,2,k*r1)*k*r1 -besselh(1,1,k*r1)*k*r1 ; ...
               -besselh(1,2,k*r2)*k*r2 -besselh(1,1,k*r2)*k*r2 ];

        b = [ 1 ; 0 ];

        R = diag(sum(abs(dM),2)); % row scaling matrix
        x = (R\dM)\(R\b);

        V = M*x;
        I = -2*pi/z;
        Za = V(1)/I;
    else
        Za = z/(2*pi)*besselh(0,2,k*r1)/(besselh(1,2,k*r1)*k*r1);
    end
        
    Z2 = [ Z2 Za ];
    
end

fhz = freqs/(2*pi);
plot( fhz, real(Z1), '-r', fhz, real(Z2), '*r', fhz, imag(Z1), '-b', fhz, imag(Z2), '*b' );
