%% Triangulation refinement
%
% One can specify how many levels of refinement desired by changing
% refine_num. Triangualtion refinement is preformed using publically
% available Matlab code by Dirk-Jan Kroon found at 
% http://www.mathworks.com/matlabcentral/fileexchange/16215-triangular-mesh-refinement/content/refinepatch.m

restoredefaultpath;
addpath(genpath('~/Desktop/cp_matrices'));
addpath('refinepatch_version2b')

%% Load original pig to be refined
[Fac, Ver] = plyread('pig_loop2.ply', 'tri');
Ver = [Ver(:,1), Ver(:,3), Ver(:,2)];  % switch axis pig is aligned on

FV.vertices = Ver;
FV.faces = Fac;

%% Refine triangulation
refine_num = 1;
for i = 1:refine_num
FV = refinepatch(FV);
end

%% Save refinement in .mat file
Faces = FV.faces;
Vertices = FV.vertices;
save(['pig_refined',num2str(refine_num),'.mat'],'Faces','Vertices');
