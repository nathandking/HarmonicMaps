%% Harmonic mapping from a cylinder to sphere.
%
% The harmonic map is computed using the closest point method as
% described in Nathan King and Steve Ruuth, "Solving variational problems 
% and partial differential equations mapping between manifolds via the 
% closest point method", Journal of Computational Physics, 2017.
%
% This example shows the computation of a harmonic map between two generic
% manifolds.

restoredefaultpath;
addpath(genpath('~/Desktop/cp_matrices'));
tic

%% Construct a grid in the embedding space
dx = 0.025;           % grid size.
% make vectors of x, y, positions of the grid
x1d = (-2.0:dx:2.0)';
y1d = x1d;
z1d = x1d;

%% Banding: do calculation in a narrow band around the manifolds
dim = 3;    % dimension
p = 3;      % interpolation degree
order = 2;  % Laplacian order
% "band" is a vector of the indices of the points in the computation
% band.  The formula for bw is found in [Ruuth & Merriman 2008] and
% the 1.0001 is a safety factor.
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));

%% Find closest points for cylinder
% meshgrid is only needed for finding the closest points, not afterwards
[xx, yy, zz] = meshgrid(x1d, y1d, z1d);
% function cpSphere for finding the closest points on a sphere
[cyl_cpx, cyl_cpy, cyl_cpz, cyl_dist] = cpCylinder(xx,yy,zz);

band = find(abs(cyl_dist) <= bw*dx);

% make into vectors
cyl_cpxg = cyl_cpx(:); cyl_cpyg = cyl_cpy(:); cyl_cpzg = cyl_cpz(:);
cyl_cpxg = cyl_cpxg(band); cyl_cpyg = cyl_cpyg(band); cyl_cpzg = cyl_cpzg(band);

%% Find closest points on the sphere
% For each point (x,y,z), store closest point on the sphere (cpx,cpy,cpz)

% meshgrid is only needed for finding the closest points, not afterwards
[xx, yy, zz] = meshgrid(x1d, y1d, z1d);
% function cpSphere for finding the closest points on a sphere
[sph_cpx, sph_cpy, sph_cpz, sph_dist] = cpSphere(xx,yy,zz);
% make into vectors
sph_cpxg = sph_cpx(:); sph_cpyg = sph_cpy(:); sph_cpzg = sph_cpz(:);
sph_cpxg = sph_cpxg(band); sph_cpyg = sph_cpyg(band); sph_cpzg = sph_cpzg(band);

%% Function u in the embedding space
% this makes u into a vector, containing only points in the band
u1_init = sph_cpxg;
u2_init = sph_cpyg;
u3_init = sph_cpzg;

%% Create Laplacian matrix for heat equation
disp('Constructing Laplacian matrix');
L = laplacian_3d_matrix(x1d,y1d,z1d, order, band, band);

%% Construct extension matrix for cylinder
disp('Constructing extension matrix for cylinder')
E_cyl = interp3_matrix(x1d, y1d, z1d, cyl_cpxg, cyl_cpyg, cyl_cpzg, p, band);

%% Load initial "texture map"
%imag = 'texture_cboard';
imag = 'four_circles';
dx_str = num2str(dx);
load(['initial_maps/',imag,'_cyl_sph_dx',dx_str(3:end),'.mat']);
C = C_sph/255;

%% Plot the initial texture mapped image.
[V, F, C_tri] = tri_color([u1_init, u2_init, u3_init], C);

% cut off triangles that close top and bottom; these triangles are only
% constructed due to the method of triangulation, these are not actually
% part of the map.
tol = 0.01;
Vz = V(:,3);
mVz = max(abs(Vz));
id_max1 = find(mVz - tol <= abs(Vz));
id_max2 = find(abs(Vz) <= mVz);
id_max = intersect(id_max1, id_max2);

Lia_max = ismember(F,id_max);
F_new = F(sum(Lia_max,2)<=1,:);

figure;
patch('Vertices', V, 'Faces', F_new,'FaceVertexCData', C_tri,'FaceColor',...
    'flat','edgecolor', 'none');
camlight(-180,120); view([-150 0]); axis off; axis equal; material dull;
set(gcf,'PaperSize',[5.25 5.25])
     
%% Add noise to map.
N1 = 0.075*randn(length(cyl_cpxg),1);
N2 = 0.075*randn(length(cyl_cpyg),1);
N3 = 0.075*randn(length(cyl_cpzg),1);
cyl_u1 = cyl_cpxg + N1;
cyl_u2 = cyl_cpyg + N2;
cyl_u3 = cyl_cpzg + N3;

[cpcyl_u1, cpcyl_u2, cpcyl_u3] = cpCylinder(cyl_u1, cyl_u2, cyl_u3);
% map noisy data back onto sphere and assign initial u1, u2, u3.
[u1, u2, u3] = cpSphere(cpcyl_u1, cpcyl_u2, cpcyl_u3);

%% Plot the noisy texture mapped initial image.
[V, F, C_tri] = tri_color([u1, u2, u3], C);

% cut off triangles that close top and bottom;
Vz = V(:,3);
mVz = max(abs(Vz));
id_max1 = find(mVz - tol <= abs(Vz));
id_max2 = find(abs(Vz) <= mVz);
id_max = intersect(id_max1, id_max2);

Lia_max = ismember(F,id_max);
F_new = F(sum(Lia_max,2)<=1,:);
    
figure;
patch('Vertices', V, 'Faces', F_new,'FaceVertexCData', C_tri,'FaceColor',...
    'flat','edgecolor', 'none');
camlight(-180,120); view([-150 0]); axis off; axis equal; material dull;
set(gcf,'PaperSize',[5.25 5.25])    
 
%% Time-stepping for harmonic mapping.
dt = 0.1*dx^2;
Tf = 300*dt;
numtimesteps = ceil(Tf/dt);
% adjust for integer number of steps
dt = Tf / numtimesteps;
tStart = tic;
for kt = 1:numtimesteps
    % explicit Euler timestepping
    utilde1 = u1 + dt*L*u1;
    utilde2 = u2 + dt*L*u2;
    utilde3 = u3 + dt*L*u3;
    
    unew1 = E_cyl*utilde1;
    unew2 = E_cyl*utilde2;
    unew3 = E_cyl*utilde3;

    [u1, u2, u3] = cpSphere(unew1, unew2, unew3);
    
    t = kt*dt
end
t_scheme = toc(tStart);

%% Plot difused mapping visualizing with the texture mapped image.
[V, F, C_tri] = tri_color([u1, u2, u3], C);

% cut off triangles that close top and bottom;
Vz = V(:,3);
mVz = max(abs(Vz));
id_max1 = find(mVz - tol <= abs(Vz));
id_max2 = find(abs(Vz) <= mVz);
id_max = intersect(id_max1, id_max2);

Lia_max = ismember(F,id_max);
F_new = F(sum(Lia_max,2)<=1,:);

figure;
patch('Vertices', V, 'Faces', F_new,'FaceVertexCData', C_tri,'FaceColor',...
    'flat','edgecolor', 'none');
camlight(-180,120); view([-150 0]); axis off; axis equal; material dull;
set(gcf,'PaperSize',[5.25 5.25])

tfin = toc;
    