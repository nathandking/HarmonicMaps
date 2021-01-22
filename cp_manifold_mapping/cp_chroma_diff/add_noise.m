%% Add color salt and pepper noise to image.
function [Cr_noise, Cg_noise, Cb_noise] = add_noise(chroma,per)

m2 = size(chroma,1);
k = ceil(per*m2);

noise_idx = randperm(m2,k);
color_idx = randi(3,k,1);

C_noise = chroma;
C_noise(noise_idx,:) = zeros(k,3);

for i = 1:k
    C_noise(noise_idx(i),color_idx(i)) = 1;
end

Cr_noise = C_noise(:,1);
Cg_noise = C_noise(:,2);
Cb_noise = C_noise(:,3);
