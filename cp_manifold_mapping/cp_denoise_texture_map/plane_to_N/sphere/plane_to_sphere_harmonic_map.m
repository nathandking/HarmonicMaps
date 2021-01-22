%% Harmonic mapping from a plane to sphere.
%
% Visualization of the mapping is shown by diffusing a noisy texture mapped
% image. The initial texture map is computed using multidimensional
% scaling described in G. Zigelman, R. Kimmel and N. Kiryati, "Texture
% mapping using surface flattening via multidimensional scaling", IEEE
% Transactions on Visualization and Computer Graphics, Vol. 8, No. 2, 2002.
%
% For further details see Nathan King and Steven Ruuth, "Solving variational 
% problems and partial differential equations mapping between manifolds via 
% the closest point method", Journal of Computational Physics, 2017.

restoredefaultpath;
addpath(genpath('~/Desktop/cp_matrices'));
tic

%% Construct a grid in the embedding space.
dx = 0.05;           % grid size.
% make vectors of x, y, positions of the grid
x1d = (-2.0:dx:2.0)';
y1d = x1d;
z1d = x1d;

%% Find closest points on the surface.
% For each point (x,y,z), store closest point on the sphere (cpx,cpy,cpz)

% meshgrid is only needed for finding the closest points, not afterwards
[xx, yy, zz] = meshgrid(x1d, y1d, z1d);
% function cpSphere for finding the closest points on a sphere
[cpx, cpy, cpz, dist] = cpSphere(xx,yy,zz);
% make into vectors
cpxg = cpx(:); cpyg = cpy(:); cpzg = cpz(:);

%% Banding: do calculation in a narrow band around the sphere.
dim = 3;    % dimension
p = 3;      % interpolation degree
order = 2;  % Laplacian order
% "band" is a vector of the indices of the points in the computation
% band.  The formula for bw is found in [Ruuth & Merriman 2008] and
% the 1.0001 is a safety factor.
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));
band = find(abs(dist) <= bw*dx);

% store closest points in the band;
cpxg = cpxg(band); cpyg = cpyg(band); cpzg = cpzg(band);

%% Take closest points as initial map value on N.
u1_init = cpxg;
u2_init = cpyg;
u3_init = cpzg;

%% Load initial mapping of image onto sphere.
init_texture_map = 'initial_maps/sphere_MDS_101.mat';
size_XS = str2double(init_texture_map(end-6:end-4));

texture_image = 'initial_maps/square_sphere_MDS_101.mat';

load(init_texture_map);
load(texture_image);
W = double(Urgb);

%% Create Laplacian matrix for heat equation
disp('Constructing Laplacian matrix');
L = laplacian_3d_matrix(x1d,y1d,z1d, order, band, band);

%% Construct an interpolation matrix for plotting.
[xp,yp,zp] = sphere(size_XS-1); % number of points matches initial map size
xp1 = xp(:); yp1 = yp(:); zp1 = zp(:);

% Eplot is a matrix which interpolations data onto the plotting grid
Eplot = interp3_matrix(x1d, y1d, z1d, xp1, yp1, zp1, p, band);

%% Plot the initial texture mapped image.
U_plot_init = [Eplot*u1_init, Eplot*u2_init, Eplot*u3_init];
[V, F, C] = tri_color(U_plot_init, W); 

figure;
patch('Vertices', V, 'Faces', F,'FaceVertexCData', C,'FaceColor',...
    'flat','edgecolor', 'none');
camlight(-90,-200); view([90 0]);
material dull;
axis([-1 1 -1 1 -1 1]); axis off; axis equal;

%% Add noise to map.
sigma = 0.075;
rng('shuffle');
N = sigma*[randn(length(u1_init),1),...
           randn(length(u2_init),1),...
           randn(length(u3_init),1)];

u1 = u1_init + N(:,1);
u2 = u2_init + N(:,2);
u3 = u3_init + N(:,3);

% map noisy data back onto sphere and assign initial u1, u2, u3
[u1, u2, u3] = cpSphere(u1, u2, u3);

%% Plot the noisy texture mapped initial image.
U_plot_noisy = [Eplot*u1, Eplot*u2, Eplot*u3];
[V, F, C] = tri_color(U_plot_noisy, W); 

figure;
patch('Vertices', V, 'Faces', F,'FaceVertexCData', C,'FaceColor',...
    'flat','edgecolor', 'none');
camlight(-90,-200); view([90 0]);
material dull;
axis([-1 1 -1 1 -1 1]); axis off; axis equal;

%% Time-stepping for harmonic mapping.
dt = 0.1*dx^2;
Tf = 60*dt;
numtimesteps = ceil(Tf/dt);
% adjust for integer number of steps
dt = Tf / numtimesteps;
tStart = tic;
for kt = 1:numtimesteps
    % explicit Euler timestepping
    unew1 = u1 + dt*L*u1;
    unew2 = u2 + dt*L*u2;
    unew3 = u3 + dt*L*u3;

    [u1, u2, u3] = cpSphere(unew1, unew2, unew3);
    
    t = kt*dt
end
t_scheme = toc(tStart);

%% plot diffused mapping visualizing with the texture mapped image.
U_plot = [Eplot*u1, Eplot*u2, Eplot*u3];
[V, F, C] = tri_color(U_plot, W); 
    
figure;
patch('Vertices', V, 'Faces', F,'FaceVertexCData', C,'FaceColor',...
    'flat','edgecolor', 'none');
camlight(-90,-200); view([90 0]);
material dull;
axis([-1 1 -1 1 -1 1]); axis off; axis equal;

tfin = toc;
