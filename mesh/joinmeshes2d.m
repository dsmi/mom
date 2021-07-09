function [ e, v ] = joinmeshes2d( e1, v1, e2, v2 )
% [ e, v ] = joinmeshes2d( e1, v1, e2, v2 )
%
% Join together two meshes
%

e = [ e1; e2 + size(v1, 1) ];
v = [ v1; v2 ];
