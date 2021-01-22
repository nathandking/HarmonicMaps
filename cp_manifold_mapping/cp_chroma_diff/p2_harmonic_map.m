%% Color image denoising via chroma diffusion using 2-harmonic map.
%
% For description of concepts see B. Tang and G. Sapiro, "Color image
% enhancement via chromaticity diffusion", IEEE Transactions on Image
% Processing, Vol. 10, No. 5, 2001.
%
% The 2-harmonic map is computed using the closest point method as
% described in Nathan King and Steve Ruuth, "Solving variational problems 
% and partial differential equations mapping between manifolds via the 
% closest point method", Journal of Computational Physics, 2017.
%
% Note: this is not the most efficient implementation possible. The
% embedding space consists of the entire cube, i.e. no banding around the
% sphere. This is done to easily obtain the planar denoised image from one
% of the z-slices.

restoredefaultpath;
addpath(genpath('~/Desktop/cp_matrices'));
tic

%% Create/load initial image.
imag = '3colors';
%imag = 'fine_fall_day';

pix = 128;
if strcmp(imag,'3colors')
    Im_orig = make_3color(pix);
else
    Im_orig = imread(['../planar_images/',imag,'.jpg']);
    Im_orig = imresize(Im_orig,[pix pix]);
    Im_orig = double(Im_orig);
end

% plot original image
figure;
imshow(uint8(Im_orig(2:pix-1,2:pix-1,:)));

%% Construct a grid in the embedding space.
dx = 4/(pix-1);           % grid size at one pixel resolution

% make vectors of x, y, positions of the grid
x1d = (-2.0:dx:2.0)';
y1d = x1d;
z1d = x1d;

%% Create hypercube extending from planar image.
[xx, yy, zz] = meshgrid(x1d, y1d, z1d);

% create band just on the outermost part of the hypercube
xg = xx(:); yg = yy(:); zg = zz(:);
Xg = abs([xg, yg, zg]);
max_row_val = max(Xg,[],2);
band = find(max_row_val <= 2-dx);

%% Calculate chroma.
Img_orig = reshape(Im_orig,[pix^2,3]);
I = sqrt(sum(abs(Img_orig).^2,2));

Cr = Img_orig(:,1)./I;
Cg = Img_orig(:,2)./I;
Cb = Img_orig(:,3)./I;

%% Add noise to chroma.
percent = 0.05;          % percent of pixels to add noise to.
[Cr_noise, Cg_noise, Cb_noise] = add_noise([Cr,Cg,Cb],percent); 
Img_noise = [I.*Cr_noise, I.*Cg_noise, I.*Cb_noise];

% plot noisy image
figure;
Im_noise = reshape(Img_noise,[pix,pix,3]);
imshow(uint8(Im_noise(2:pix-1,2:pix-1,:)));

%% Replicate chroma onto all closest for every constant z-plane.
C1 = Cr_noise(:,ones(size(zz,3),1));
C2 = Cg_noise(:,ones(size(zz,3),1));
C3 = Cb_noise(:,ones(size(zz,3),1));

u1 = C1(band);
u2 = C2(band);
u3 = C3(band);

%% Create Laplacian matrix for heat equation
disp('Constructing Laplacian matrix')
L = laplacian_3d_matrix(x1d,y1d,z1d, 2, band);

%% Time evolution for 2-harmonic mapping.
dt = 0.1*dx^2;        
numtimesteps = 30;
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

%% Find chroma of just one z-slice of hypercube.
c1 = ones(pix^3,1);
c2 = ones(pix^3,1);
c3 = ones(pix^3,1);

c1(band) = u1;
c2(band) = u2;
c3(band) = u3;

portion = pix^2+1:2*pix^2;
Img_denoised = [I.*c1(portion), I.*c2(portion), I.*c3(portion)];

% plot denoised image
figure;
Im_denoised = reshape(Img_denoised,[pix,pix,3]);
imshow(uint8(Im_denoised(2:pix-1,2:pix-1,:)));

tfin = toc;