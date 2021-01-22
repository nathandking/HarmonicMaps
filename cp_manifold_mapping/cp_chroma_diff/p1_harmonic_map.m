%% Color image denoising via chroma diffusion using 1-harmonic map.
%
% For description of concepts see B. Tang and G. Sapiro, "Color image
% enhancement via chromaticity diffusion", IEEE Transactions on Image
% Processing, Vol. 10, No. 5, 2001.
%
% The 1-harmonic map is computed using the closest point method as
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

%% Create derivative matrices.
disp('Constructing upwind scheme matrices for first derivatives')
[Dxb,Dxf,Dyb,Dyf,Dzb,Dzf] = firstderiv_upw1_3d_matrices(x1d,y1d,z1d,band);

%% Time evolution for 1-harmonic mapping.
dt = 0.5*dx^2;
numtimesteps = 125;

norm_Jac_u = zeros(length(u1),1);
delta = 1e-16;
tStart = tic;
for kt = 1:numtimesteps
    grad_u1 = [Dxf*u1, Dyf*u1, Dzf*u1];
    grad_u2 = [Dxf*u2, Dyf*u2, Dzf*u2];
    grad_u3 = [Dxf*u3, Dyf*u3, Dzf*u3];
 
    Jac_u_vec = [grad_u1, grad_u2, grad_u3];
    norm_Jac_u = sqrt(sum(Jac_u_vec.^2,2)) + delta;
    
    div1 = Dxb*(grad_u1(:,1)./norm_Jac_u) + Dyb*(grad_u1(:,2)./norm_Jac_u)...
           + Dzb*(grad_u1(:,3)./norm_Jac_u);
    div2 = Dxb*(grad_u2(:,1)./norm_Jac_u) + Dyb*(grad_u2(:,2)./norm_Jac_u)...
           + Dzb*(grad_u2(:,3)./norm_Jac_u);
    div3 = Dxb*(grad_u3(:,1)./norm_Jac_u) + Dyb*(grad_u3(:,2)./norm_Jac_u)...
           + Dzb*(grad_u3(:,3)./norm_Jac_u);       
    
    % forward Euler timestepping
    unew1 = u1 + dt*div1;
    unew2 = u2 + dt*div2;
    unew3 = u3 + dt*div3;

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