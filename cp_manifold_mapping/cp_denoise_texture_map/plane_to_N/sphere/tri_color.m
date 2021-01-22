%% Triangulation of surface using Delaunay triangulation.
function [TriP, Tri_n, Tri_Urgb] = tri_color(XS, Urgb)

% determine the unique set of points.
[uni_XS, ia, ~] = unique(XS,'rows'); 
uni_Urgb = Urgb(ia,:);

DT = delaunayTriangulation(uni_XS(:,1),uni_XS(:,2),uni_XS(:,3));
[Tri, Tri_points] = freeBoundary(DT);

% map Tri_points onto XS.
[~,ib,ic]=intersect(uni_XS, Tri_points, 'rows');

% extract the needed colors using the ib map.
Tri_Urgb = uni_Urgb(ib,:);

% permute the surface triangulation points using ic map
TriP = Tri_points(ic,:);

% map the point numbers used in triangle definitions
% NOTE: for that inverse map is needed
iic(ic) = 1:length(ic);
Tri_n = iic(Tri);

