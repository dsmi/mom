%
% "Overclocked capacitor" - impedance of the waveguide formed by a pair of
% round plates. Compare the results against the analytical solution.
%

addpath(genpath([ pwd, '/..' ]));

% The geometry - round with a hole in center which is the port.
inch2meter = 2.54e-2;
r2 = 0.9*inch2meter; % 0.9 inches radius
r1 = r2*3e-2;
d = 10e-3*inch2meter; % 10 mils separation
n1 = 8;
n2 = 100;
eps1 = eps0*4.3;

[ e1, v1 ] = mkcir2d(r1, n1);
e1(:, [ 1 2 ]) = e1(:, [ 2 1 ]);

[ e2, v2 ] = mkcir2d(r2, n2);
e = [ e1; e2 + size(v1, 1) ];
v = [ v1; v2 ];

% The port edges
c1 = find_edges2d(e, v, 0, 0, r1*1.1);
ports = { c1' };

%% plotmesh2d(e,v,ports,1);

Ya = []; % Analytical solution
Y1 = []; % Simulated
Y2 = []; % Simulated via admittance

fhz0 = 1e7;
fhz1 = 1e10;
%freqs = logspace(log10(fhz0*2*pi),log10(fhz1*2*pi),20);
freqs = linspace(fhz0*2*pi,fhz1*2*pi,500);

%freqs = 850e6*2*pi; % first resonance

for freq = freqs,

    freq

    % Parameters of the plane.
    Yplane = j*freq*eps1/d;
    Zplane = j*freq*mu0*d;
    k = sqrt(-Yplane*Zplane);

    wavelen = 2*pi/k

    fintgsl = @(rsrc,robs) intg_helmsl2d(k,rsrc,robs);
    fintgdl = @(rsrc,robs) intg_helmdl2d(k,rsrc,robs);

    [ Z, u, q ] = extractz2(e, v, ports, fintgsl, fintgdl);
    Z = Z*Zplane;
    Y1 = cat(3, Y1, inv(Z));

    Y=extracty2(e, v, ports, fintgsl, fintgdl)/Zplane;
    Z=1/Y;
    Y2 = cat(3, Y2, inv(Z));

    % Analytical solution
    M = [ besselj(0, k*r1) bessely(0, k*r1) ; ...
		                            besselj(0, k*r2) bessely(0, k*r2) ];
    %M = [ besselh(0,2,k*r1) besselh(0,1,k*r1) ; ...
    %      besselh(0,2,k*r2) besselh(0,1,k*r2) ]

    dM = [ -besselj(1, k*r1)*k*r1 -bessely(1, k*r1)*k*r1 ; ...
		                                         -besselj(1, k*r2)*k*r2 -bessely(1, k*r2)*k*r2 ];
    %dM = [ -besselh(1,2,k*r1)*k -besselh(1,1,k*r1)*k ; ...
    %       -besselh(1,2,k*r2)*k -besselh(1,1,k*r2)*k ]

    b = [ 1 ; 0 ];
    x = dM\b;

    V = M*x;
    I = -2*pi/Zplane;
    Z = V(1)/I;
    Ya = cat(3, Ya, inv(Z));

end

%% %% img = mkimage(e, v, fintgsl, fintgdl, u, q, 128, 128);

% Capacitor approximation
C = eps1*pi*r2*r2/d;
Yc = j*freqs*C;

% Estimate inductance
%L = 1.5e-10;
L = 9.5682e-011;%L = 1.7813e-010%1.0160e-010
Yl = 1./(j*freqs*L);
Ylc = 1./(1./Yc + 1./Yl); % series connection

% Ladder network
nr = 1000;
dr = (r2 - r1)/nr;
r = linspace(r1 + dr/2, r2 - dr/2, nr);
Cr = eps1*pi*((r + dr/2).^2 - (r - dr/2).^2)/d;
Lr = mu0*d./(2*pi*r)*dr;
Yladder = 0;
for i=nr:-1:1
    Yladder = Yladder + j*freqs*Cr(i);
    Zladder = 1./Yladder + j*freqs*Lr(i);
    Yladder = 1./Zladder;
end

fhz = freqs/(2*pi);
semilogy(fhz,abs(shiftdim(Ya)),'-xg', fhz,abs(Yc),'-xb',fhz,abs(Ylc),'-or',fhz,abs(Yladder),'-vm')
xlabel('freq, Hz');
ylabel('Im(Y), S');
legend('Ya', 'Yc', 'Ylc', 'Yladder');

%% %% k = freqs.*sqrt(mu0*eps1);



tswrite('round2.y1p', freqs/(2*pi), Y1)
tswrite('round2a.y1p', freqs/(2*pi), Ya)
tswrite('round2lc.y1p', freqs/(2*pi), shiftdim(Ylc, -1))

Za = 1./Ya;
Ltrace = 0.666e-9*2;
Zt = Za + shiftdim(j*freqs*Ltrace, -1);
Yt = 1./Zt;
tswrite('round2trace.y1p', freqs/(2*pi), Yt)

Yc = j*freqs*C;
tswrite('round2c.y1p', freqs/(2*pi), shiftdim(Yc, -1))
