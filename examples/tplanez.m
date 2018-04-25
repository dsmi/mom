%
% Transmission plane with two ports.
%

addpath(genpath([ pwd, '/..' ]));

% The geometry
inch2meter = 2.54e-2;
d = 30*1e-3*inch2meter; % 30 mils separation
n1 = 10;
rport = 0.005*inch2meter;
nr = 10;
maxl = 0.2*inch2meter; % max edge len
% Dielectric params
lt = 0.02;
er0 = 4.3;
fr = 1e9;
% metal conductivity
sigma = 5.8e7;

v0 = [ 0 0 ; 0 2 ; 2.5 2 ; 2.5 1 ; 3 1 ; 3 0 ; 0 0 ]*inch2meter;
%v0 = [ 0 0 ; 0 2 ; 3 2 ; 3 0 ; 0 0 ]*inch2meter;
e0 = [ (1:size(v0,1)-1)' (2:size(v0,1))' ];

% Split edges (if needed)
[ e, v ] = subdiv2d(e0, v0, maxl);

% Remove duplicated vertices, which may be created when splitting edges
[ e, v ] = rmdups2d(e, v);

% First port
[ e1, v1 ] = mkcir2d(rport, nr);
v1 = v1 + repmat([ 1 1 ]*inch2meter, size(v1, 1), 1);
e = [ e; e1 + size(v, 1) ];
v = [ v; v1 ];

% Second port
[ e2, v2 ] = mkcir2d(rport, nr);
v2 = v2 + repmat([ 2 1 ]*inch2meter, size(v2, 1), 1);
e = [ e; e2 + size(v, 1) ];
v = [ v; v2 ];

e(:, [ 1 2 ]) = e(:, [ 2 1 ]); % flip inside out to get correct normals

% The port edges
c1 = find_edges2d(e, v, 1*inch2meter, 1*inch2meter, rport*1.1);
c2 = find_edges2d(e, v, 2*inch2meter, 1*inch2meter, rport*1.1);
ports = { c1' c2' };

%plotmesh2d(e,v,ports,1);

Z1 = []; % Simulated impedance

fhz0 = 1e8;
fhz1 = 5e9;
freqs = logspace(log10(fhz0*2*pi),log10(fhz1*2*pi),20);
%freqs = linspace(fhz0*2*pi,fhz1*2*pi,1000);

for freq = freqs,
    % Parameters of the plane.
    er = debye(er0, lt, fr, freq/(2*pi));
    Yplane = j*freq*eps0*er/d;
    Zs = 2*sqrt(j*freq*mu0/sigma); % surface impedance
    Zplane = Zs + j*freq*mu0*d;
    k = sqrt(-Yplane*Zplane);

    fintgsl = @(rsrc,robs) intg_helmsl2d(k,rsrc,robs);
    fintgdl = @(rsrc,robs) intg_helmdl2d(k,rsrc,robs);

    [ Z, u, q ] = extractz2(e, v, ports, fintgsl, fintgdl);
    Z = Z*Zplane;
    Z1 = cat(3, Z1, Z);
end

img = mkimage(e, v, fintgsl, fintgdl, u(:,1), q(:,1), 60, 40);
imgx=linspace(0,3,7);
imgy=linspace(2,0,5);
imagesc(imgx, imgy, imag(img*Zplane), [ -10 10 ])
colorbar
%% xlabel('$x$, inch');
%% ylabel('$y$, inch');
%% set (gcf,'paperposition',[ 0 0 5.5 3.6 ])
%% print('-dtex', 'vplane.tex');

%% hldata = dlmread('tplane2.z2p', ' ');
%% fhl = hldata(:, 1);
%% Z11 = hldata(:, 2) + j*hldata(:, 3);
%% Z12 = hldata(:, 4) + j*hldata(:, 5);

%% fhz = freqs/(2*pi);
%% plot(fhz,imag(shiftdim(Z1(1,2,:))),'-r','LineWidth',2,...
%%      fhz,imag(shiftdim(Z1(1,1,:))),'-m','LineWidth',2,...
%%      fhl,imag(Z12),'-b','LineWidth',2,...
%%      fhl,imag(Z11),'-c','LineWidth',2)
%% grid on
%% axis([0,fhz1,-60,60])
%% xlabel('$freq, \hertz$');
%% ylabel('$Im(Z), \ohm$');
%% legend('$Z_{1,2}$', '$Z_{1,1}$', '$Z_{1,2}$-HL', '$Z_{1,1}$-HL');
%% ylim([-35 35])
%% NumTicks = 6;
%% L = get(gca,'XLim');
%% set(gca,'XTick',linspace(L(1),L(2),NumTicks))
%% set (gcf,'paperposition',[ 0 0 5.5 3.6 ])
%% print('-dtex', 'zplane.tex');
