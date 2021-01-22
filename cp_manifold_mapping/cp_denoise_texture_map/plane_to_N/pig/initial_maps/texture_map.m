%% Construction of texture map
% This code uses MDS flattened coordinates to map planar image on the 
% surface of an object.

%% Load flattened coordinates
tri = 'pig_refined1';
load([tri,'_MDS.mat']);
newx = newxy(:,1); newy = newxy(:,2);

%% Load texture image
imag = 'geoArt';

T = imread(['../../../../planar_images/',imag,'.jpg']);
T = flipud(T);
T = rot90(T,-1);

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

[uu,vv] = ndgrid(u,v);

%% Interpolate colors in (u,v) coordinates to (newx,newy) coordinates
Ur = griddata(uu(:), vv(:), Tr, newxy(:,1), newxy(:,2),'nearest');
Ug = griddata(uu(:), vv(:), Tg, newxy(:,1), newxy(:,2),'nearest');
Ub = griddata(uu(:), vv(:), Tb, newxy(:,1), newxy(:,2),'nearest');

Ur = Ur/max(Ur);
Ug = Ug/max(Ug);
Ub = Ub/max(Ub);

UUr = 233*ones(n,1)/255;
UUg = 150*ones(n,1)/255;
UUb = 122*ones(n,1)/255;

Urgb = [UUr, UUg, UUb];
Urgb(sub_idx,:) = [Ur, Ug, Ub];

tfin_color = toc;
save([imag,'_',tri,'_texture_map.mat'],'Urgb','tfin_color');