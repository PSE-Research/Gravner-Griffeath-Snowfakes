function ppm2png(ppm_path)
I = imread(ppm_path);
ppm_size = size(I);
ppm_L = ppm_size(1);

% 30deg shear
shear = [...
1	0.5 	0;...
0   1	0; ...
0   0	1
];
tform = affine2d(shear);
% fill background with color in point (2,2,1)
J = imwarp(I, tform, 'FillValues', I(2,2,1));

y_add = tan(pi/180*30) * ppm_L;
png = imcrop(J, [0, y_add/2, ppm_L, ppm_L]);
imshow(png)
