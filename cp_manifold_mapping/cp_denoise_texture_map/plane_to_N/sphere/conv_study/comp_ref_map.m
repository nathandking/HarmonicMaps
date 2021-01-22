%% Computation of reference map for plane to sphere convergence study
restoredefaultpath;
addpath(genpath('~/Desktop/cp_matrices'));
tic

%% Construct a grid in the embedding space
dx = 0.025; % grid spacing

max_c = 1.0 + 3*dx*bw; 
max_N = ceil(max_c/dx);
max_cube = dx*max_N;
% make vectors of x, y, positions of the grid
x1d = (-max_cube:dx:max_cube)';
y1d = x1d;
z1d = x1d;

%% Find closest points on the surface
% For each point (x,y,z), store closest point on the surface (cpx,cpy,cpz)

% meshgrid is only needed for finding the closest points, not afterwards
[xx, yy, zz] = meshgrid(x1d, y1d, z1d);
% function cpSphere for finding the closest points on a sphere
[cpx, cpy, cpz, dist] = cpSphere(xx,yy,zz);
% make into vectors
cpxg = cpx(:); cpyg = cpy(:); cpzg = cpz(:);
clear cpx cpy cpz;

%% Compute bandwidth
dim = 3;    % dimension
p = 3;      % interpolation degree
order = 2;  % Laplacian order
% "band" is a vector of the indices of the points in the computation
% band.  The formula for bw is found in [Ruuth & Merriman 2008] and
% the 1.0001 is a safety factor.
bw = 1.0001*sqrt((dim-1)*((p+1)/2)^2 + ((order/2+(p+1)/2)^2));

%% Banding: do calculation in a narrow band around the sphere
band = find(abs(dist) <= bw*dx);
% store closest points in the band;
cpxg = cpxg(band); cpyg = cpyg(band); cpzg = cpzg(band);

%% Take initial map to be closest points; no noise added.
u1 = cpxg; 
u2 = cpyg; 
u3 = cpzg;

%% Create Laplacian matrix for heat equation
disp('Constructing Laplacian matrix');
L = laplacian_3d_matrix(x1d,y1d,z1d, order, band, band);

%% Construct an interpolation matrix for to put points on texture map locations
% plotting grid on sphere, based on parameterization
size_XS = 401;
[xp,yp,zp] = sphere(size_XS-1);   % number of points must match initial map size.
xp1 = xp(:); yp1 = yp(:); zp1 = zp(:);

disp('Constructing interpolation matrix')
Eplot = interp3_matrix(x1d, y1d, z1d, xp1, yp1, zp1, p, band);
    
%% Time-stepping for harmonic mapping.
dt = 0.1*dx^2;
Tf = 0.015;
numtimesteps = ceil(Tf/dt);
% adjust for integer number of steps
dt = Tf / numtimesteps;
tStart = tic;
for kt = 1:numtimesteps
    % explicit Euler timestepping
    u1 = u1 + dt*L*u1;
    u2 = u2 + dt*L*u2;
    u3 = u3 + dt*L*u3;

    [u1, u2, u3] = cpSphere(u1, u2, u3);
    
    t = kt*dt
end
t_scheme = toc(tStart);

%% harmonic map closest points on texture map locations.
U_ref = [Eplot*u1, Eplot*u2, Eplot*u3];
U = [u1, u2, u3];
dx_str = num2str(dx);
Tf_str = num2str(Tf);
tfin = toc;

save(['ref_map_dx',dx_str(3:end),'_Tf',Tf_str(3:end),'_p',num2str(p),'.mat'],'U_ref','tfin');

