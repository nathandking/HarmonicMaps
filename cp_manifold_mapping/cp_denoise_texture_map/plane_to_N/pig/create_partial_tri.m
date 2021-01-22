%% Compute partial triangulation from an entire triangulated object
%
% Returns a .txt file with partial triangulation vertices and faces.

%% Load entire triangulated object
tri_mesh = 'pig_refined1';
load([tri_mesh,'.mat']);

xp = Vertices(:,1);
yp = Vertices(:,2);   % switch plotting axis
zp = Vertices(:,3);

Xp = [xp, yp, zp];

%% Plot entire triangulation
figure(1); clf;
patch('Vertices', Vertices, 'Faces', Faces,'FaceColor', 'none');
xlabel('x')
ylabel('y')
zlabel('z')
view([110,10])

%% Index partial triangulation to obtain vertices
x_g = find(xp >= 1/10);
z_g = find(zp >= 1/8);
y_g = find(yp >= -1/2); 

y_l = find(yp <= 1/2);
z_l = find(zp <= 1);
 
x_intersect = intersect(x_g, z_g);
z_intersect = intersect(y_l, z_l);

xz_intersect = intersect(x_intersect, z_intersect);
sub_idx = intersect(y_g,xz_intersect);

XS_subset = Xp(sub_idx,:);

%% Map Tri_points onto XS to determine corresponding faces
[V,ia,ib] = intersect(XS_subset, Xp, 'rows');

Lia = ismember(Faces,ib);
Liaa = ismember(Lia,[1,1,1],'rows');
Tri = Faces(Liaa == 1,:);

iic(ib) = 1:length(ib);
F = iic(Tri);

%% Plot partial mesh
figure(2); clf;
patch('Vertices', V, 'Faces', F,'FaceColor', 'none');
xlabel('x')
ylabel('y')
zlabel('z')
view([110,10])

%% Write vertices and faces of partial triangulation to .txt file
dlmwrite(['partial_',tri_mesh,'.txt'],[size(V,1), size(F,1)],'delimiter','\t','precision',10);
dlmwrite(['partial_',tri_mesh,'.txt'],V,'-append','delimiter','\t','precision',10);
dlmwrite(['partial_',tri_mesh,'.txt'],F-ones(size(F)),'-append','delimiter','\t');
