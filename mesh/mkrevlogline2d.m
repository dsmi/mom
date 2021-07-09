function [ e, v ] = mkrevlogline2d(dx, nx)
% [ e, v ] = mkrevlogline2d(dx, nx)
%
% Makes a right-to-left logarithmically spaced line parallel
% to x axis, dimensions are (-dx/2 - dx/2)
%

v = [ flip( -( (logspace(0, 1, nx+1)-1)*dx/9-dx/2 )', 1 ) repmat(0, nx+1, 1) ];
e = [ (1:nx)' (2:nx+1)' ];
