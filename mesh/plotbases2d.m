function plotbases2d(edges,verts,bases)

vx = verts(:,1);
vy = verts(:,2);
plot(vx(edges)', vy(edges)', 'b-');
xlabel('x');
ylabel('y');

r1 = verts(edges(:,1),:);
r2 = verts(edges(:,2),:);
ev = r2 - r1;

l = sqrt(sum(ev.^2,2));
t = ev ./ l(:,ones(1,2));
n = [ t(:,2) -t(:,1) ];

% Bases scale (auto-calculate)
s = max( max(verts) - min(verts) ) * 1e-2;

hold on

inedges = bases( :, 1 );

b = r1(inedges,:); % beginning
e = r2(inedges,:); % end
p = e+n(inedges,:)*s; % pinnacle
plot( transpose( [ b(:,1) p(:,1) ] ), transpose( [ b(:,2) p(:,2) ] ), 'r-');
plot( transpose( [ p(:,1) e(:,1) ] ), transpose( [ p(:,2) e(:,2) ] ), 'r-');

outedges = bases( :, 2 );

b = r1(outedges,:); % beginning
e = r2(outedges,:); % end
p = b+n(outedges,:)*s; % pinnacle
plot( transpose( [ b(:,1) p(:,1) ] ), transpose( [ b(:,2) p(:,2) ] ), 'g-');
plot( transpose( [ p(:,1) e(:,1) ] ), transpose( [ p(:,2) e(:,2) ] ), 'g-');


hold off

