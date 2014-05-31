function [ e, v ] = mkline2d(dx, nx)
% [ e, v ] = mkline2d(dx, nx)
%
% Makes a line parallel to x axis, dimensions are (-dx/2 - dx/2)
%

v = [ linspace(-dx/2, dx/2, nx+1)' repmat(0, nx+1, 1) ];
e = [ (1:nx)' (2:nx+1)' ];
