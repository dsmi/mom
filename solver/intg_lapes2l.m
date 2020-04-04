function v = intg_lapes2l( edgelen, epsout, epsin, srce, ua, ub, obse, uc, ud )
% v = intg_lapes2l( edgelen, epsout, epsin, srce, ua, ub, obse, uc, ud )
%
% Calculates self-terms (where the source and observation edges coincide)
% of the E matrix which transforms the boundary element charges to
% discontinuity of the normal component of displacement.
%
% Where f(u) and f(v) are the weighting functions, r is the
% source-to-observation vector, R = |r|.
% Depending on the need/context, it ca be viewed as:
%  (a) Electric field due to the charge on the source edge, dotted with
%      the observation normal (can be used for enforcing continuity of
%      the electric field at the dielectric-dielectric interface)
%  (b) Potential at the source due to the double electric layer at
%      observation.
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

% Do the calculation for all the edges -- for simplicity
eout = epsout(srce);
ein  = epsin(srce);
edlen = edgelen(srce);
diffe = eout-ein;
diffe( find( diffe == 0 ) ) = 1.0; % Fix for the conductor edges
intff = 1/6*edlen.*(ua.*(2*uc+ud)+ub.*(2*ud+uc));
vall = (1/(2*eps0))*intff.*(eout+ein)./diffe;
                          
% Find singular pairs
sidx = find( srce == obse );

% And then save the results for singular edges only
v(sidx) = vall(sidx);                      

%% % Integral of the normal derivative is zero for the self terms -- due
%% % to the symmetry.
%% v(sidx) = 0*v(sidx);
