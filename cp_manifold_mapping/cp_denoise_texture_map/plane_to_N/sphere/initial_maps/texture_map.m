%% Construction of texture map
% This code uses MDS flattened coordinates to map planar image on the 
% surface of an object.

%% Load flattened coordinates
load('sphere_MDS_101.mat');
n = 100;
newx = newxy(:,1); newy = newxy(:,2);

%% Load texture image
image = 'square';
T = imread(['../../../../planar_images/',image,'.jpg']);

T = flipud(T);
Tr = double(T(:,:,1));
Tg = double(T(:,:,2));
Tb = double(T(:,:,3));

Tr = Tr(:);
Tg = Tg(:);
Tb = Tb(:);

%% Scale texture image to lie within (newx,newy) coordinates
xmin = min(newx); ymin = min(newy);
xmax = max(newx); ymax = max(newy);

u = (xmin:(xmax-xmin)/(size(T,1)-1):xmax)';
v = (ymin:(ymax-ymin)/(size(T,2)-1):ymax)';

%% Interpolate colors in (u,v) coordinates to (newx,newy) coordinates
[uu,vv] = ndgrid(u,v);
Ur = griddata(uu(:), vv(:), Tr, newxy(:,1), newxy(:,2),'natural');
Ug = griddata(uu(:), vv(:), Tg, newxy(:,1), newxy(:,2),'natural');
Ub = griddata(uu(:), vv(:), Tb, newxy(:,1), newxy(:,2),'natural');

Ur = Ur/max(Ur);
Ug = Ug/max(Ug);
Ub = Ub/max(Ub);

Urgb = [Ur, Ug, Ub];

save([image,'_sphere_MDS_',num2str(n+1),'.mat'],'Urgb');
