%% "Texture map" for mapping between cylinder and sphere.
% The psuedo texture map between the cylinder and sphere is created by
% mapping points on the sphere to points on the cylinder by taking the 
% spheres radial direction.
%
% An image is place on the cylinder simply by wrapping the planar image
% into a cylinder. The image can then be put on the sphere using the map
% between the cylinder and sphere discussed above.

restoredefaultpath;
addpath(genpath('~/Desktop/cp_matrices'));

%% Set up embedding space
dx = 0.025;  %grid size

% make vectors of x, y, positions of the grid
x1d = (-2.0:dx:2.0)';
y1d = x1d;
z1d = x1d;

%% Banding: do calculation in a narrow band around the sphere
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

%% Find for each closest point on unit sphere its corresponding point 
% radially outward on the unit radius cylinder

z_cyl = sph_cpzg./sqrt(1-sph_cpzg.^2);
x_cyl = sqrt(1+z_cyl.^2).*sph_cpxg;
y_cyl = sqrt(1+z_cyl.^2).*sph_cpyg;

%% Flattening of the cylinder to a plane
[x_pix,~] = cart2pol(x_cyl,y_cyl);
y_pix = z_cyl;

%% Pixel coordinates of a flattened cylinder
imag = 'four_circles';
Img = imread(['../../../planar_images/',imag,'.jpg']);
pix = size(Img,1);
Img = Img(1:pix,1:pix,:);

Im = 256*ones(3*pix, 3*pix, 3);
Im(2*pix+1:3*pix,pix+1:2*pix,:) = Img;

C = double(reshape(Im,(3*pix)^2,3));
Cr = C(:,1);
Cg = C(:,2);
Cb = C(:,3);

R = 1;
th = linspace(-pi,pi,3*pix);
z = linspace(-2.0,2.0,3*pix);

[x_pixel, y_pixel] = ndgrid(th,z);

%% Interpolate color values on image cylinder to sphere_cylinder
Fr = scatteredInterpolant(x_pixel(:), y_pixel(:), Cr, 'linear');
Fg = scatteredInterpolant(x_pixel(:), y_pixel(:), Cg, 'linear');
Fb = scatteredInterpolant(x_pixel(:), y_pixel(:), Cb, 'linear');

% evaluate scattered interpolant function.
Cr_sph = Fr(x_pix, y_pix);
Cg_sph = Fg(x_pix, y_pix);
Cb_sph = Fb(x_pix, y_pix);

C_sph = [Cr_sph, Cg_sph, Cb_sph];

M = [x_cyl, y_cyl, z_cyl];
N = [sph_cpxg, sph_cpyg, sph_cpzg];
dx_str = num2str(dx);
save([imag,'_cyl_sph_dx',dx_str(3:end),'.mat'],'M','N','C_sph')
 