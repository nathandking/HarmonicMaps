%% Code to make the three color (R,G,B) image
function [Im]=make_3color(pix)

Im = zeros(pix, pix, 3);
% assign left half of image to be red (maybe add some variation after)
Im(:,1:pix/2,1) = 256*ones(pix,pix/2);
% assign right upper quarter to be green
Im(1:pix/2,pix/2+1:pix,2) = 256*ones(pix/2,pix/2);
% assign right upper quarter to be green
Im(pix/2+1:pix,pix/2+1:pix,3) = 256*ones(pix/2,pix/2);