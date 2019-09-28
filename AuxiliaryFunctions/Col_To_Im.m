function mI = Col_To_Im(mP, vImage_size, vPatch_size)

patch_height = vPatch_size(1);              
patch_width  = vPatch_size(2);

image_height = vImage_size(1);
image_width  = vImage_size(2);

mI = zeros(vImage_size);

k = 1;
for x = 1 : (image_width - patch_width) + 1
    for y = 1 : (image_height - patch_height) + 1
        mI(y:y+patch_height-1, x:x+patch_width-1) = ...
         mI(y:y+patch_height-1, x:x+patch_width-1) + reshape(mP(:,k), vPatch_size);
        k = k + 1;
    end
end

M = ones(vImage_size - vPatch_size + 1);
H = ones(vPatch_size);

Mask = conv2(M, H);

mI = mI ./ Mask;