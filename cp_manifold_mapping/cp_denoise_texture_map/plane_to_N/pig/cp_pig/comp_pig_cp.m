%% Computation of the closest point function for triangulated pig
restoredefaultpath;
addpath(genpath('~/Desktop/cp_matrices'));

%% Load mat file contains the triangles
disp('loading pig');
tri = 'pig_refined1';
load(['../',tri,'.mat']);
%% Computation of closest points uses mex implemention for efficiency
%mex tri2cp_helper.c

dx = 0.05;   % grid size
disp('running tri2cp');
[IJK,DIST,CP,XYZ,WHICH_FACES] = tri2cp(Faces, Vertices, dx, -2);

%% Save closest point function in .mat file
dx_str = num2str(dx);
save([tri,'_CP_dx',dx_str(3:end),'.mat'],'IJK','DIST','CP','XYZ');
