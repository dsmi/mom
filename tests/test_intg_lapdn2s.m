function test_intg_lapdn2s
% 
%% addpath(genpath([ pwd, '/..' ]));

% Simple case -- edge len = 0.5, both basis functions = 1 everywhere
v = intg_lapdn2s( 0.5, 2.0, 3.0, 1, 1, 1, 1, 1, 1 );
vt = 1/(2*eps0) * 0.5 * (2.0+3.0)/(2.0-3.0);
assertEquals(vt, v, 1e-3) % value is ~1e11

% edge len = 0.5, both basis functions = change from 0 to 1
v = intg_lapdn2s( 0.5, 2.0, 3.0, 1, 0.0, 1, 1, 0.0, 1 );
vt = 1/(2*eps0) * 1/6 * (2.0+3.0)/(2.0-3.0);
assertEquals(vt, v, 1e-3) % value is ~1e11

% edge len = 2, basis functions 0 to 1 and 1 to 0
v = intg_lapdn2s( 2.0, 2.0, 3.0, 1, 0.0, 1, 1, 1, 0.0 );
vt = 1/(2*eps0) * 1/3 * (2.0+3.0)/(2.0-3.0);
assertEquals(vt, v, 1e-3) % value is ~1e11
