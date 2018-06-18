function S1 = renorms(S0, Z0, Z1)
% S1 = renorms(S0, Z0, Z1)
%
% Renormalize S-parameters
%

R = diag((Z1-Z0)./(Z1+Z0));
A = diag(sqrt(Z1./Z0)*1./(Z1+Z0));

S1 = S0*0; % pre-allocate

for i=1:size(S0,3)
  S = S0(:,:,i);
  size(inv(A)*(S-R)*inv(eye(size(S)) - R*S)*A);
  S1(:,:,i) = inv(A)*(S-R)*inv(eye(size(S)) - R*S)*A;
end
