[img, info] = mha_read_volume('./036.mha');
lbl = mha_read_volume('./ours_036_c1.mha');

% downsampling
factor = 2;
img = img(1:factor:end, 1:factor:end, 1:factor:end);
lbl = lbl(1:factor:end, 1:factor:end, 1:factor:end);
dims = size(img);

% data format
img = double(img);
img = 254 * (img - min(img(:))) / (max(img(:)) - min(img(:)));
lbl = 255 * (lbl >= 0.7*65535);
img = uint8(img);
lbl = uint8(lbl);
img = img(:, :, end:-1:1);
lbl = lbl(:, :, end:-1:1);

% decompose into list of cubes
cube_size = 48;
ita = 0;
fold = dims / cube_size + ita;
ovlap = ceil((fold * cube_size - dims) ./ (fold - 1));
fold = ceil(((fold - 1) .* ovlap + dims) ./ cube_size);

% get indexes
r_s_list = [];
r_e_list = [];
c_s_list = [];
c_e_list = [];
h_s_list = [];
h_e_list = [];
for H = 0:fold(3)-1
    h_s = H*cube_size - H*ovlap(3) + 1;
    h_e = h_s + cube_size - 1;
    if h_e > dims(3)
        h_s = dims(3) - cube_size + 1;
        h_e = h_s + cube_size - 1;
    end
    for C = 0:fold(2)-1
        c_s = C*cube_size - C*ovlap(2) + 1;
        c_e = c_s + cube_size - 1;
        if c_e > dims(2)
            c_s = dims(2) - cube_size + 1;
            c_e = c_s + cube_size - 1;
        end
        for R = 0:fold(1)-1
            r_s = R*cube_size - R*ovlap(1) + 1;
            r_e = r_s + cube_size - 1;
            if r_e > dims(1)
                r_s = dims(1) - cube_size + 1;
                r_e = r_s + cube_size - 1;
            end
            r_s_list = [r_s_list, r_s];
            r_e_list = [r_e_list, r_e];
            c_s_list = [c_s_list, c_s];
            c_e_list = [c_e_list, c_e];
            h_s_list = [h_s_list, h_s];
            h_e_list = [h_e_list, h_e];
        end
    end
end

% convolutional animation
info.Dimensions = dims;
info.PixelDimensions = info.PixelDimensions * factor;
info.DataType = 'uchar';
img_cov = img;
mha_write(img_cov(:,:,end:-1:1), info, './1000_animation.mha');
for ii = 1:length(r_s_list)
    r_s = r_s_list(ii);
    r_e = r_e_list(ii);
    c_s = c_s_list(ii);
    c_e = c_e_list(ii);
    h_s = h_s_list(ii);
    h_e = h_e_list(ii);
    
    img_cov(r_s:r_e, c_s:c_e, h_s:h_e) = 0.5 * img_cov(r_s:r_e, c_s:c_e, h_s:h_e);
    
    img_cov(r_s, c_s, h_s:h_e) = 51;
    img_cov(r_s, c_e, h_s:h_e) = 51;
    img_cov(r_e, c_s, h_s:h_e) = 51;
    img_cov(r_e, c_e, h_s:h_e) = 51;
    
    img_cov(r_s, c_s:c_e, h_s) = 51;
    img_cov(r_s, c_s:c_e, h_e) = 51;
    img_cov(r_e, c_s:c_e, h_s) = 51;
    img_cov(r_e, c_s:c_e, h_e) = 51;
    
    img_cov(r_s:r_e, c_s, h_s) = 51;
    img_cov(r_s:r_e, c_s, h_e) = 51;
    img_cov(r_s:r_e, c_e, h_s) = 51;
    img_cov(r_s:r_e, c_e, h_e) = 51;
    
    mha_write(img_cov(:,:,end:-1:1), info, ['./data/', num2str(1000+ii), '_animation.mha']);
    
    img_cov(r_s:r_e, c_s:c_e, h_s:h_e) = img_cov(r_s:r_e, c_s:c_e, h_s:h_e) ...
                                         .* (lbl(r_s:r_e, c_s:c_e, h_s:h_e) / 255);
%     img_cov(r_s:r_e, c_s:c_e, h_s:h_e) = lbl(r_s:r_e, c_s:c_e, h_s:h_e);
%     img_cov(r_s:r_e, c_s:c_e, h_s:h_e) = zeros([cube_size, cube_size, cube_size]);
end
