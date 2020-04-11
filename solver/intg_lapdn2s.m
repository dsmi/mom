function v = intg_lapdn2s( edgelen, epsout, epsin, srce, ua, ub, obse, uc, ud )
% v = intg_lapdn2s( edgelen, epsout, epsin, srce, ua, ub, obse, uc, ud )
%
% Calculates self-terms (where the source and observation edges coincide)
% of the E matrix which transforms the boundary element charges to
% discontinuity of the normal component of displacement.
%
%  The self-terms are given by:
%    e(i,j) = \intg 1/(2*eps0) * (epso+epsi)/(epso-epsi) f_i(r) * f_j(r) dr
%  where epso and epsi are the dielectric permittivities inside and outside.
%
%  edgelen       - length of the edges
%  epsout, epsin - dielectric permittivity outside and inside of an edge
%  srce, ua, ub  - source edge index and f(u) at the endpoints
%  obse, ua, ub  - observation edge index and f(u) at the endpoints
% 
%  v - the resulting self-terms, N-by-1
%

% Number of edges
N = size(srce,1);

% This pre-allocates the output array
v = 0 * ua;

% Find singular pairs
sidx = find( srce == obse );

% Do the calculation for the singular edges only
eout = epsout( srce(sidx) );
ein  = epsin( srce(sidx) );
elen = edgelen( srce(sidx) );
a = ua( sidx );
b = ub( sidx );
c = uc( sidx );
d = ud( sidx );

% Integral of the product of the basis functions
intff = 1/6*elen.*(a.*(2*c+d)+b.*(2*d+c));

% And then save the results for singular edges
v(sidx) = (1/(2*eps0))*intff.*(eout+ein)./(eout-ein+1e-300);
