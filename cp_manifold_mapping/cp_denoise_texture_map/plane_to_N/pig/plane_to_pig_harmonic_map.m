%% Harmonic mapping from a plane to pig
%
% Visualization of the mapping is shown by diffusing a noisy texture mapped
% image. The initial texture map is computed using multidimensional
% scaling described in G. Zigelman, R. Kimmel and N. Kiryati, "Texture
% mapping using surface flattening via multidimensional scaling", IEEE
% Transactions on Visualization and Computer Graphics, Vol. 8, No. 2, 2002.
%
% The closest point function is computed from the triangulated pig using
% code in cp_matrices by S. Ruuth, C. Macdonald and/or students. The 
% required code can be found in ./cp_pig.
%
% For further details of this example see Nathan King, "Solving variational
% problems and partial differential equations mapping between manifolds via
% the closest point method", M.Sc. Thesis, Simon Fraser University, 2015.

restoredefaultpath;
addpath(genpath('~/Desktop/cp_matrices'));
tic

%% Load triangulated pig and texture map
tri = 'pig_refined1';
init_map = [tri,'.mat'];
load(init_map);

imag = 'geoArt';
load(['initial_maps/',imag,'_',tri,'_texture_map.mat']);
C = double(Urgb);

%% Load cp function for pig and compute band
dx = 0.05;
str_dx = num2str(dx);
cp_file = ['cp_pig/',tri,'_CP_dx',str_dx(3:end),'.mat'];
load(cp_file);

i = IJK(:,1);
j = IJK(:,2);
k = IJK(:,3);
dist = DIST;
cpxg = CP(:,1);
cpyg = CP(:,2);
cpzg = CP(:,3);
xg = XYZ(:,1);
yg = XYZ(:,2);
zg = XYZ(:,3);

x1d=-2.0:dx:2.0;
y1d=x1d;
z1d=x1d;
nx=length(x1d);
ny=length(y1d);
nz=length(z1d);

dim = 3;
p = 3;  % degree interp
order = 2;  % laplacian order, griddata hardcoded for 2

%%sanity checks
xtest = x1d(1) + (i-1)*dx;
ytest = y1d(1) + (j-1)*dx;
ztest = z1d(1) + (k-1)*dx;
if ( (norm(xtest - xg, inf) > 1e-14) ||  ...
        (norm(ytest - yg, inf) > 1e-14) || ...
        (norm(ztest - zg, inf) > 1e-14) )
    error('sanity fail');
end

% here is one place where meshgrid comes in: note ordering here
band = sub2ind([ny,nx,nz], j,i,k);

%% Construct an interpolation matrix for plotting
xp1 = Vertices(:,1); 
yp1 = Vertices(:,2); 
zp1 = Vertices(:,3);

% Eplot is a matrix which interpolations data onto the plotting grid
disp('Constructing plotting interpolation matrix');
Eplot = interp3_matrix(x1d, y1d, z1d, xp1, yp1, zp1, p, band);
    
%% Create Laplacian matrix for heat equation
disp('Constructing Laplacian matrix');
L = laplacian_3d_matrix(x1d,y1d,z1d, order, band, band);

%% Plot original texture mapped image on plotting points  
Original = [Eplot*cpxg, Eplot*cpyg, Eplot*cpzg];

figure;
patch('Vertices', Original, 'Faces', Faces,'FaceVertexCData', C,'FaceColor',...
    'flat','EdgeColor', 'none');
camlight(55,15); view([65,-10]); axis off; axis([-0.9 0.9 -0.9 0.9 0 0.8])
material dull; set(gcf,'PaperSize',[4.5 4]);

%% Add normally distributed noise to closest points
sigma = 0.075;
N = sigma*[randn(size(CP,1),1), randn(size(CP,1),1), randn(size(CP,1),1)];
U_noise = CP + N;

%% Try to do CP_N(U_noise)
% Best idea is to do knnsearch onto plotting vertices
[IDX] = knnsearch(Vertices,U_noise,'NSMethod','exhaustive');

% initial mapped points on computational domain
u1 = Vertices(IDX,1);
u2 = Vertices(IDX,2);
u3 = Vertices(IDX,3);

%% Plot noisy map on plotting vertices
noisy_texture_map_plot(Vertices, Faces, C, sigma);

%% Construct scattered interpolant function
Fx = scatteredInterpolant(xg, yg, zg, cpxg, 'linear');
Fy = scatteredInterpolant(xg, yg, zg, cpyg, 'linear');
Fz = scatteredInterpolant(xg, yg, zg, cpzg, 'linear');

%% Time-stepping for harmonic map from plane to pig
dt = 0.1*dx^2;
Tf = 17*dt;
numtimesteps = ceil(Tf/dt);
% adjust for integer number of steps
dt = Tf / numtimesteps;
tStart = tic;
for kt = 1:numtimesteps
    % explicit Euler timestepping
    unew1 = u1 + dt*L*u1;
    unew2 = u2 + dt*L*u2;
    unew3 = u3 + dt*L*u3;

    % evaluate scattered interpolant function to determine closest points
    u1 = Fx(unew1, unew2, unew3);
    u2 = Fy(unew1, unew2, unew3);
    u3 = Fz(unew1, unew2, unew3);

    t = kt*dt
end
t_scheme = toc(tStart);

%% Plot diffused mapping visualizing with the texture mapped image
U_plot = [Eplot*u1, Eplot*u2, Eplot*u3];

figure;
patch('Vertices', U_plot, 'Faces', Faces,'FaceVertexCData', C,'FaceColor',...
    'flat','EdgeColor', 'none');
camlight(55,15); view([65,-10]); axis off; axis([-0.9 0.9 -0.9 0.9 0 0.8])
material dull; set(gcf,'PaperSize',[4.5 4]);

tfin = toc;
