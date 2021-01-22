%% Multidimensional scaling used to flatten an object into the plane.
%
% Method is described in G. Zigelman, R. Kimmel and N. Kiryati, "Texture
% mapping using surface flattening via multidimensional scaling", IEEE
% Transactions on Visualization and Computer Graphics, Vol. 8, No. 2, 2002.
%
% The pairwise geodesic distances are computed using .c code in the
% geo_dist/ folder. This publically available code is by Danil Kirsanov 
% found at https://code.google.com/p/geodesic/. It is not the most 
% efficient method for geodesic distance computation.
%
% Note: this code performs the MDS on a subset of triangulated pig.

%% Load triangulated pig vertices to be flattened via MDS
tri = 'pig_refined1';
load(['../',tri,'.mat']);

n = size(Vertices,1);
xp = Vertices(:,1);
yp = Vertices(:,2);
zp = Vertices(:,3);
Xp = [xp, yp, zp];

%% Determine subset of pig to be flattened
x_g = find(xp >= 1/10);
z_g = find(zp >= 1/8);
y_g = find(yp >= -1/2); 

y_l = find(yp <= 1/2);
z_l = find(zp <= 1);
 

x_intersect = intersect(x_g, z_g);
z_intersect = intersect(y_l, z_l);

xz_intersect = intersect(x_intersect, z_intersect);
ssub_idx = intersect(y_g,xz_intersect);

XS_subset = Xp(ssub_idx,:);

%% Map Tri_points onto XS.
[V,ia,ib] = intersect(XS_subset, Xp, 'rows');
sub_idx = ssub_idx(ia);

%% Load geodesic distance between pairs of points
M = load(['partial_',tri,'_geo_dist.txt']);
M = 0.5*(M+M');  % average the doublely computed distances

%% Apply classical multidimensional scaling to flatten coordinates
J = eye(size(M,1)) - ones(size(M,1))./(size(M,1));
B = -0.5 * J * M * J;
B = 0.5*(B+B');       % In theory B should be symmetric, so fix this.
clear J;

% Find largest eigenvalues + their eigenvectors
opts.tol = 1e-8;
[Q, L,~] = eigs(B,2,'LA',opts);
clear B;

% Extract the coordinates.
newxy = [sqrt(L(1,1)).*Q(:,1), sqrt(L(2,2)).*Q(:,2)];
clear L Q;

save([tri,'_MDS.mat'],'newxy','sub_idx','n')
