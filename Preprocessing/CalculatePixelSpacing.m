%% Determine the pixel spacing from the microscope slide thickness in view
% The No. 1 slides are 0.13 to 0.17 mm in thickness:
% https://us.vwr.com/store/product/4645817/vwr-micro-cover-glasses-rectangular

close all

% Slide thickness (mm)
slide_thickness_mm = mean([.13 .17]);

% Load a BG image
bg_image_file = 'D:\GitHub\FlyTripod_eLife_2021\Preprocessing\180325SideView_BG.jpg';
bg_image = imread(bg_image_file);
bg_image = rgb2gray(bg_image);
bg_image = medfilt2(bg_image, [20, 20]);
bg_image = edge(bg_image,'sobel');

f = figure;
imshow(bg_image)

% Rectangle drawn around the microscope slide [x y w h]:
rectangle_dims = [0.0717863105175293 0.693606755126658 0.853088480801335 0.0832328106151982];
annotation(f,'rectangle',...
    rectangle_dims,...
    'Color',[1 0 0]);

% The height is the slide thickness in pixels
slide_thickness_pixels = rectangle_dims(4);

% Determine the pixel spacing (pixels/mm)
pixels_mm = slide_thickness_pixels / slide_thickness_mm;
