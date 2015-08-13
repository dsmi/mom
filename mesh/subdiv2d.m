function [ e, v ] = subdiv2d(e, v, maxl)
% [ e, v ] = subdiv2d(e, v, maxl)
%
% Subdivides/splits edges longer than maxl
%

% Edge beginnings and ends
r1 = v(e(:,1), :);
r2 = v(e(:,2), :);

% Calculate edge lengths.
el = sqrt(sum((r2 - r1).^2,2)); % Column vector of length N.

% Split the edges if too long. Vectorize this?
newe = [];
for ie=1:size(e,1),
	n = ceil(el(ie)/maxl-1e-6);
	if n > 1,
		x = linspace(r1(ie,1), r2(ie,1), n+1);
		y = linspace(r1(ie,2), r2(ie,2), n+1);
		newe = [ newe; size(v,1) + [ (1:n)' (2:n+1)' ] ];
		v = [ v; [ x' y' ] ];
	else
		newe = [ newe; e(ie,:) ];
	endif
end

e = newe;
