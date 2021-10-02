%
% Here we calculate two-port Y-parameters of the ring-shaped cavity section
% with symmetric voltage and current. We then calculate Z-parameters of the
% "overclocked capacitor" as in roundz by connecting two ring sections -- outer
% one which is terminated and inner two-port.
%

addpath(genpath([ pwd, '/..' ]));

% One-port Z-parameters as in roundz.m
function Z = rz1( Yplane, Zplane, r1, r2 )
    
    k = sqrt(-Yplane*Zplane);
    
    % Analytical solution
    M = [ besselj(0, k*r1) bessely(0, k*r1) ; ...
		                            besselj(0, k*r2) bessely(0, k*r2) ];

    dM = [ -besselj(1, k*r1)*k*r1 -bessely(1, k*r1)*k*r1 ; ...
		                                         -besselj(1, k*r2)*k*r2 -bessely(1, k*r2)*k*r2 ];

    b = [ 1 ; 0 ];
    x = dM\b;

    V = M*x;
    I = -2*pi/Zplane;
    Z = V(1)/I;
    
end    

% Two-port Y-parameters
function Y = ry( Yplane, Zplane, r1, r2 )
    
    k = sqrt(-Yplane*Zplane);
    
    % Analytical solution
    M = [ besselj(0, k*r1) bessely(0, k*r1) ; ...
          besselj(0, k*r2) bessely(0, k*r2) ];

    dM = [ -besselj(1, k*r1)*k*r1 -bessely(1, k*r1)*k*r1 ; ...
           -besselj(1, k*r2)*k*r2 -bessely(1, k*r2)*k*r2 ];

    v = [ 1 0 ; 0 1 ];
    b = M\v;
    Y = [ -1 0 ; 0 1 ]*dM*b*2*pi/Zplane;

end    

% The geometry - round with a hole in center which is the port.
inch2meter = 2.54e-2;
r3 = 0.9*inch2meter; % 0.9 inches radius
r2 = r3*3e-1;
r1 = r3*3e-2;
d = 10e-3*inch2meter; % 10 mils separation
eps1 = eps0*4.3;
Yc = 1e4;
Zc = 1e-4;

Ya = [ ]; % Analytical solution
Yb = [ ]; % Analytical solution of two parts

fhz0 = 1e7;
fhz1 = 1e10;
freqs = logspace(log10(fhz0*2*pi),log10(fhz1*2*pi),20);
%% freqs = 1e3*2*pi

for freq = freqs,

    freq

    % Parameters of the plane.
    Yplane = Yc+j*freq*eps1/d;
    Zplane = j*freq*mu0*d;
    k = sqrt(-Yplane*Zplane);

    wavelen = 2*pi/k

    Z = rz1( Yplane, Zplane, r1, r3 );
    Ya = [ Ya inv(Z) ];

    % Inner two-port section.
    Y = ry( Yplane, Zplane, r1, r2 );

    % Outer 'terminated' section
    Z2 = inv( ry( Yplane, Zplane, r2, r3 ) );
    Yl = inv( Z2(1,1) );
    
    V1 = (Y(2,2)+Yl)/((Y(2,2)+Yl)*Y(1,1)-Y(2,1)*Y(1,2));
    Yb = [ Yb 1./V1 ];
    
end

plot( freqs, imag(Ya), '-r', freqs, imag(Yb), 'b*' )


%% freq = 1e10

%% % Dielectric params
%% lt = 0.02;
%% er0 = 4.3;
%% fr = 1e9;

%% % metal conductivity
%% sigma = 5.8e7;

%% % Parameters of the plane.
%% er = debye( er0, lt, fr, freq/(2*pi) );
%% Yplane = j*freq*eps0*er/d;
%% Zs = 2*sqrt(j*freq*mu0/sigma); % surface impedance
%% Zplane = Zs + j*freq*mu0*d;
%% k = sqrt(-Yplane*Zplane)

%% Y = ry( Yplane, Zplane, r1, r2 )

