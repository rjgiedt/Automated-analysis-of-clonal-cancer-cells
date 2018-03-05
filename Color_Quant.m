% Program for analyzing color intensities.  Files to analyze should be
% placed in 3 seperate folders according to color (red, blue, green).
% the directory paths for each color should then be input into the "Load
% Files block for each color respectively.  Each file should have the same
% number of images, and the image order should each correspond to each
% other.  

%% Setup 
clear all; clc; close all;

%% User defined values
diff_value = 1000;

%% Load files 
red_images = ReadImgs('/Users/randygeidt/Desktop/red','*.tif');
green_images = ReadImgs('/Users/randygeidt/Desktop/green','*.tif');
blue_images = ReadImgs('/Users/randygeidt/Desktop/blue','*.tif');

red_length = length(red_images(:,512,512));
green_length = length(green_images(:,512,512));
blue_length = length(blue_images(:,512,512));

if red_length ~= green_length && green_length ~= blue_length
error('Must have equal number of images in each folder')
end

number_images = red_length;

%% Prompt for filtered video output name
output_video_fname = 'Output.avi';
[output_video_fname,output_video_path] = uiputfile('*.avi','Input output video file name');

if isequal(output_video_fname,0) || isequal(output_video_path,0)
   disp('User selected Cancel')
   return;
else
   fprintf('User selected "%s" as output video file\n',output_video_fname);
end

%% Construct a VideoWriter object and view its properties. 
writerObj = VideoWriter(output_video_fname);

%Set the frame rate to 1 frames per second.
frame_rate = 1;
writerObj.FrameRate = frame_rate;  
open(writerObj);

blank = zeros(512,512);
  blank(:,:) = NaN;


%% Analyze images
for n = 1:1%number_images
n
    
initial_red = red_images(n,:,:);
initial_green = green_images(n,:,:);
initial_blue = blue_images(n,:,:);

red = squeeze(initial_red);
green = squeeze(initial_green);
blue = squeeze(initial_blue);

figure; imagesc(red); axis square; caxis([0 1000]);
figure; imagesc(green); axis square; caxis([0 1000]);
figure; imagesc(blue); axis square; caxis([0 1000]);

red_original = red;
green_original = green;
blue_original = blue;

red = uint16(red);
green = uint16(green);
blue = uint16(blue);

%% Find which images have color for use in segmentation
% Calculate image maximum intensities
[value_red, location_red] = max(red_original(:));
[value_green, location_green] = max(green_original(:));
[value_blue, location_blue] = max(blue_original(:));

% Calculate average intensity for each color image
red_average = mean(red_original(:));
green_average = mean(green_original(:));
blue_average = mean(blue_original(:));

% Calculate intensity differences 
red_difference = value_red - red_average;
green_difference = value_green - green_average;
blue_difference = value_blue - blue_average;

% Red - Green - Blue
if red_difference > diff_value && green_difference > diff_value &&...
        blue_difference > diff_value
    
    rgbImage_pre = imfuse(imadjust(red), imadjust(green), 'blend');
    rgbImage = imfuse(rgbImage_pre, imadjust(blue), 'blend');
% rgbImage(:,:,1) = imadjust(red);
% rgbImage(:,:,2) = imadjust(green);
% rgbImage(:,:,3) = imadjust(blue);

% Red
elseif red_difference > diff_value && green_difference < diff_value &&...
        blue_difference < diff_value
    rgbImage = imadjust(red);
    
% Green
elseif red_difference < diff_value && green_difference > diff_value &&...
        blue_difference < diff_value
    
    rgbImage = imadjust(green);
   
% Blue    
elseif red_difference < diff_value && green_difference < diff_value &&...
        blue_difference > diff_value

    rgbImage = imadjust(blue);
  
    
% Red-Green
elseif red_difference > diff_value && green_difference > diff_value &&...
        blue_difference < diff_value
 rgbImage = imfuse(imadjust(red), imadjust(green), 'blend');
%     rgbImage(:,:,1) = imadjust(red);
%     rgbImage(:,:,3) = imadjust(green);
rgbImage = uint16(rgbImage);
    
% Blue-Green
elseif red_difference < diff_value && green_difference > diff_value &&...
        blue_difference > diff_value
  rgbImage = imfuse(imadjust(green), imadjust(blue), 'blend');
%     rgbImage(:,:,1) = imadjust(green);
%     rgbImage(:,:,3) = imadjust(blue);
%     rgbImage(:,:,3) = blank;
rgbImage = uint16(rgbImage);
% Red-Blue
elseif red_difference > diff_value && green_difference < diff_value &&...
        blue_difference > diff_value
    
    rgbImage = imfuse(imadjust(red), imadjust(blue), 'blend');
%     rgbImage(:,:,1) = imadjust(red);
%     rgbImage(:,:,3) = imadjust(blue);
%     rgbImage(:,:,3) = blank;
    rgbImage = uint16(rgbImage);
% Nothing works
else 
    disp('Loop Skipped - No color match');
    n
    continue
end

%% Create Gray RGB
test_size = size(size(rgbImage));

if test_size(2) > 2
    rgbImage2 = rgb2gray(rgbImage);  
else
   rgbImage2 = rgbImage; 
end

%% Convert into thresholding image
test_size = size(size(rgbImage));

if test_size(2) > 2
    thresholding_image = im2double(rgb2gray(rgbImage));  
else
   thresholding_image = im2double(rgbImage); 
end

%% Segment Cells/ Find ROIs
grayBand_reverse = imadjust(thresholding_image, [0 1], [1,0]);
threshold_image = MinimaxAT(grayBand_reverse);

%% Filter Image
se = strel('disk',1); % Set disk size
B = threshold_image;
Io = imopen(B,se);
Ie = imerode(B, se);
Iobr = imreconstruct(Ie, B);
Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
Iobrcbr = double(Iobrcbr);
Iobrcbr_corrected = imadjust(Iobrcbr, [0 1], [1 0]);
Iobrcbr_filled = imfill(Iobrcbr_corrected, 'holes');

%% Post-process filter based on size and shape properties of regions
min_cutoff = 50; % Original 75... prof 45
max_cutoff = 2250; %Original 1000...prof 1500
size_filtering = bwlabel(Iobrcbr_filled);
s = regionprops(size_filtering, 'All');

%Filter based on total area sizes
area_values = [s.Area];
idx = find((min_cutoff <= area_values) & (max_cutoff >= area_values));
final_image2 = ismember(size_filtering, idx);

%Filter based on object solidity
min_cutoff_solidity = .1;
solidity_filtering = bwlabel(final_image2);
solidity_values = [s.Solidity];
idx2 = find((min_cutoff_solidity <= solidity_values));
final_image3 = ismember(final_image2, idx2);
final_image4 = double(final_image3);

% Calculate rough background values for use in filtering
background_region_rough = imadjust(final_image4, [0 1], [1 0]);
background_label_rough = bwlabel(background_region_rough);
background_intensity_rough = regionprops(background_label_rough,...
    rgbImage2, 'MeanIntensity');
background_value_rough(n) = background_intensity_rough(1,1).MeanIntensity;

% Filter based on region intensity compared to background
intensity_cutoff = 0; %Original 100
intensity_filtering = bwlabel(final_image3);
t = regionprops(intensity_filtering, rgbImage2, 'MeanIntensity');
intensity_values = [t.MeanIntensity];
corrected_intensity_values = intensity_values-background_value_rough(1,1);
idx3 = find(intensity_cutoff <= corrected_intensity_values);
final_image5 = ismember(intensity_filtering, idx3);

%Override Intensity Filtering
final_image5 = final_image4;

%% Label Areas 
threshold_image_final2 = im2double(final_image5);
[L num(n)] = bwlabel(final_image5);

%% Find RGB Values for each identified area
stats_red = regionprops(L, red_original, 'MeanIntensity');
stats_green = regionprops(L, green_original, 'MeanIntensity');
stats_blue = regionprops(L, blue_original, 'MeanIntensity');

% Find each color's background image for rough analysis
red_background = regionprops(background_label_rough,...
    red_original, 'MeanIntensity');
green_background = regionprops(background_label_rough,...
    green_original, 'MeanIntensity');
blue_background = regionprops(background_label_rough,...
    blue_original, 'MeanIntensity');

% Find each color's background smoothed image and create new intensity
% images

red_smoothed_background_15 = imopen(red_original,strel('disk',15));
green_smoothed_background_15 = imopen(green_original,strel('disk',15));
blue_smoothed_background_15 = imopen(blue_original,strel('disk',15));

red_smoothed_background_30 = imopen(red_original,strel('disk',70));
green_smoothed_background_30 = imopen(green_original,strel('disk',70));
blue_smoothed_background_30 = imopen(blue_original,strel('disk',70));

red_subtracted_15 = red_original - red_smoothed_background_15;
green_subtracted_15 = green_original - green_smoothed_background_15;
blue_subtracted_15 = blue_original - blue_smoothed_background_15;

red_subtracted_30 = red_original - red_smoothed_background_30;
green_subtracted_30 = green_original - green_smoothed_background_30;
blue_subtracted_30 = blue_original - blue_smoothed_background_30;

% background_value_rough_red(n) = mean(red_smoothed_background_15(:));
% background_value_rough_green(n) = mean(green_smoothed_background_15(:));
% background_value_rough_blue(n) = mean(blue_smoothed_background_15(:));

background_value_rough_red(n) = red_background(1,1).MeanIntensity;
background_value_rough_green(n) = green_background(1,1).MeanIntensity;
background_value_rough_blue(n) = blue_background(1,1).MeanIntensity;

% Find background subtracted intensities
stats_red_b15 = regionprops(L, red_subtracted_15, 'MeanIntensity');
stats_green_b15 = regionprops(L, green_subtracted_15, 'MeanIntensity');
stats_blue_b15 = regionprops(L, blue_subtracted_15, 'MeanIntensity');

stats_red_b30 = regionprops(L, red_subtracted_30, 'MeanIntensity');
stats_green_b30 = regionprops(L, green_subtracted_30, 'MeanIntensity');
stats_blue_b30 = regionprops(L, blue_subtracted_30, 'MeanIntensity');

%% Write Video to record regions
% Create overlay of borders
rgb_overlay = imadjust(rgbImage2);
perimeter = bwperim(L);
overlay1 = imoverlay(rgb_overlay, perimeter, [.3 1 .3]);

% Write Cell Borders Output Video    
 inputframe = im2double(overlay1);
 frame = inputframe;
 writeVideo(writerObj,frame); 

%% Compile all region intensity data
    for i=1:num(n)
    stats_red_list{n,i} = (stats_red(i).MeanIntensity);
    stats_green_list{n,i} = (stats_green(i).MeanIntensity);
    stats_blue_list{n,i} = (stats_blue(i).MeanIntensity);
    
    stats_red_list_b15{n,i} = (stats_red_b15.MeanIntensity);
    stats_green_list_b15{n,i} = (stats_green_b15.MeanIntensity);
    stats_blue_list_b15{n,i} = (stats_blue_b15.MeanIntensity);
    
    stats_red_list_b30{n,i} = (stats_red_b30.MeanIntensity);
    stats_green_list_b30{n,i} = (stats_green_b30.MeanIntensity);
    stats_blue_list_b30{n,i} = (stats_blue_b30.MeanIntensity);
    
    % Add all intensities for each for use with normalization
    summation_list{n,i} = (stats_red_list{n,i}+stats_green_list{n,i}+...
        stats_blue_list{n,i});
    summation_list_b15{n,i} = (stats_red_list_b15{n,i}+...
        stats_green_list_b15{n,i}+stats_blue_list_b15{n,i});
    summation_list_b30{n,i} = (stats_red_list_b30{n,i}+...
        stats_green_list_b30{n,i}+stats_blue_list_b30{n,i});
    
    % Create normalized values
    red_normalized{n,i} = stats_red_list{n,i}/summation_list{n,i};
    green_normalized{n,i} = stats_green_list{n,i}/summation_list{n,i};
    blue_normalized{n,i} = stats_blue_list{n,i}/summation_list{n,i};
    
    % Create smoothed background analyzed values
    red_normalized_b15{n,i} = stats_red_list_b15{n,i}/summation_list_b15{n,i};
    green_normalized_b15{n,i} = stats_green_list_b15{n,i}/summation_list_b15{n,i};
    blue_normalized_b15{n,i} = stats_blue_list_b15{n,i}/summation_list_b15{n,i};
    
    red_normalized_b30{n,i} = stats_red_list_b30{n,i}/summation_list_b30{n,i};
    green_normalized_b30{n,i} = stats_green_list_b30{n,i}/summation_list_b30{n,i};
    blue_normalized_b30{n,i} = stats_blue_list_b30{n,i}/summation_list_b30{n,i};
    
    % Create rough background analyzed values 
    stats_background_red{n,i} = (stats_red(i).MeanIntensity-background_value_rough_red(n));
    stats_background_green{n,i} = (stats_green(i).MeanIntensity-background_value_rough_green(n));
    stats_background_blue{n,i} = (stats_blue(i).MeanIntensity-background_value_rough_blue(n));
    
    % Create rough normalized background corrected values
    background_summation{n,i} = (stats_background_red{n,i}+stats_background_green{n,i}+...
        stats_background_blue{n,i});
    background_normalized_red{n,i} = stats_background_red{n,i}/background_summation{n,i};
    background_normalized_green{n,i} = stats_background_green{n,i}/background_summation{n,i};
    background_normalized_blue{n,i} = stats_background_blue{n,i}/background_summation{n,i};
  
 % Calculate hues for no background subtraction
    % red-yellow hue
    if red_normalized{n,i} >= green_normalized{n,i} &&...
            green_normalized{n,i} >= blue_normalized{n,i}
        
        hue{n,i} = 60*(green_normalized{n,i}-blue_normalized{n,i})/...
            (red_normalized{n,i}-blue_normalized{n,i});    
        
        background_normalized_hue{n,i} = 60*(background_normalized_green{n,i}-background_normalized_blue{n,i})/...
            (background_normalized_red{n,i}-background_normalized_blue{n,i});  
        
    % yellow-green hue
    elseif green_normalized{n,i} > red_normalized{n,i} &&...
            red_normalized{n,i} >= blue_normalized{n,i}
    
        hue{n,i} = 60*(2-(red_normalized{n,i}-blue_normalized{n,i})/...
            (green_normalized{n,i}-blue_normalized{n,i}));   
        
        background_normalized_hue{n,i} = 60*(2-(background_normalized_red{n,i}-background_normalized_blue{n,i})/...
            (background_normalized_green{n,i}-background_normalized_blue{n,i})); 
        
    % Green-Cyan hue    
    elseif green_normalized{n,i} >= blue_normalized{n,i} &&...
            blue_normalized{n,i} > red_normalized{n,i}
        
        hue{n,i} = 60*(2+(blue_normalized{n,i}-red_normalized{n,i})/...
            (green_normalized{n,i}-red_normalized{n,i}));   
        
        background_normalized_hue{n,i} = 60*(2+(background_normalized_blue{n,i}-background_normalized_red{n,i})/...
            (background_normalized_green{n,i}-background_normalized_red{n,i}));
        
    % Cyan-Blue hue        
    elseif blue_normalized{n,i} > green_normalized{n,i} &&...
            green_normalized{n,i} > red_normalized{n,i}
        
      hue{n,i} = 60*(4-(green_normalized{n,i}-red_normalized{n,i})/...
            (blue_normalized{n,i}-red_normalized{n,i}));      
        
      background_normalized_hue{n,i} = 60*(4-(background_normalized_green{n,i}-background_normalized_red{n,i})/...
            (background_normalized_blue{n,i}-background_normalized_red{n,i}));   
    
        % Blue-Magenta hue    
    elseif blue_normalized{n,i} > red_normalized{n,i} &&...
            red_normalized{n,i} >= green_normalized{n,i}
        
        hue{n,i} = 60*(4+(red_normalized{n,i}-green_normalized{n,i})/...
            (blue_normalized{n,i}-green_normalized{n,i}));  
        
        background_normalized_hue{n,i} = 60*(4+(background_normalized_red{n,i}-background_normalized_green{n,i})/...
            (background_normalized_blue{n,i}-background_normalized_green{n,i}));
 
    % Magenta-Red hue    
    elseif red_normalized{n,i} >= blue_normalized{n,i} &&...
            blue_normalized{n,i} > green_normalized{n,i}
        
       hue{n,i} = 60*(6-(blue_normalized{n,i}-green_normalized{n,i})/...
            (red_normalized{n,i}-green_normalized{n,i}));  
        
       background_normalized_hue{n,i} = 60*(6-(background_normalized_blue{n,i}-background_normalized_green{n,i})/...
            (background_normalized_red{n,i}-background_normalized_green{n,i}));
        
    end
    end
    
    
% Clear Data for next loop
rgbImage = [];

close all 

end   
    
%% End writer object
close(writerObj);

%% Create single variables for both hue and RGB

% Create raw values variable
rgb_raw(:,1) = cell2mat(stats_red_list(:));
rgb_raw(:,2) = cell2mat(stats_green_list(:));
rgb_raw(:,3) = cell2mat(stats_blue_list(:));

rgb_raw_b15(:,1) = cell2mat(stats_red_list_b15(:));
rgb_raw_b15(:,2) = cell2mat(stats_green_list_b15(:));
rgb_raw_b15(:,3) = cell2mat(stats_blue_list_b15(:));

rgb_raw_b30(:,1) = cell2mat(stats_red_list_b30(:));
rgb_raw_b30(:,2) = cell2mat(stats_green_list_b30(:));
rgb_raw_b30(:,3) = cell2mat(stats_blue_list_b30(:));

% Create normalized values variable
rgb_normalized(:,1) = cell2mat(red_normalized(:));
rgb_normalized(:,2) = cell2mat(green_normalized(:));
rgb_normalized(:,3) = cell2mat(blue_normalized(:));

rgb_normalized_b15(:,1) = cell2mat(red_normalized_b15(:));
rgb_normalized_b15(:,2) = cell2mat(green_normalized_b15(:));
rgb_normalized_b15(:,3) = cell2mat(blue_normalized_b15(:));

rgb_normalized_b30(:,1) = cell2mat(red_normalized_b30(:));
rgb_normalized_b30(:,2) = cell2mat(green_normalized_b30(:));
rgb_normalized_b30(:,3) = cell2mat(blue_normalized_b30(:));

% Create background raw variable
background_corrected_raw(:,1) = cell2mat(stats_background_red(:));
background_corrected_raw(:,2) = cell2mat(stats_background_green(:));
background_corrected_raw(:,3) = cell2mat(stats_background_blue(:));

% Create background normalized variable
background_normalized(:,1) = cell2mat(background_normalized_red(:));
background_normalized(:,2) = cell2mat(background_normalized_green(:));
background_normalized(:,3) = cell2mat(background_normalized_blue(:));

% Create Hues variable
Hue(:,1) = cell2mat(hue(:));
Background_Hue = cell2mat(background_normalized_hue(:));


%% Post-hoc Intensity Filtering for none-background normalized data

red_min_intensity = 100;
green_min_intensity = 800;
blue_min_intensity = 530;

for n = 1:length(rgb_raw)
if rgb_raw(n,1) <= red_min_intensity || rgb_raw(n,2) <= ...
        green_min_intensity || rgb_raw(n,3) <= blue_min_intensity 
    
   rgb_raw_filtered(n,1) = NaN;
   rgb_raw_filtered(n,2) = NaN;
   rgb_raw_filtered(n,3) = NaN;
   
   rgb_normalized_filtered(n,1) = NaN;
   rgb_normalized_filtered(n,2) = NaN;
   rgb_normalized_filtered(n,3) = NaN;
   
   Hue_filtered(n,1) = NaN;
else
rgb_raw_filtered(n,1) = rgb_raw(n,1);
rgb_raw_filtered(n,2) = rgb_raw(n,2);  
rgb_raw_filtered(n,3) = rgb_raw(n,3);

rgb_normalized_filtered(n,1) = rgb_normalized(n,1);
rgb_normalized_filtered(n,2) = rgb_normalized(n,2);
rgb_normalized_filtered(n,3) = rgb_normalized(n,3);

Hue_filtered(n,1) = Hue(n,1);
end
end

%% Post-hoc background intensity analysis
% Create Intensities variable

for n = 1:length(rgb_raw)
if background_corrected_raw(n,1) <= 0 
    filter_values_raw(n,1) = 0;
else
    filter_values_raw(n,1) = background_corrected_raw(n,1);
end

if background_corrected_raw(n,2) <= 0 
    filter_values_raw(n,2) = 0;
else
    filter_values_raw(n,2) = background_corrected_raw(n,2);
end

if background_corrected_raw(n,3) <= 0 
    filter_values_raw(n,3) = 0;
else
    filter_values_raw(n,3) = background_corrected_raw(n,3);
end
end

% Calculate intensities
for n=1:length(filter_values_raw)
intensity(n) = (filter_values_raw(n,1) + filter_values_raw(n,2) + ...
    filter_values_raw(n,3))/3;
end

% Eliminate sections not meeting minimum intensity and/or empty sets
min_filtering_intensity = 0;

for n = 1:length(rgb_raw)
if intensity(n) <= min_filtering_intensity
    final_filtered(n,1) = NaN;
    final_filtered(n,2) = NaN; 
    final_filtered(n,3) = NaN;
else
    final_filtered(n,1) = filter_values_raw(n,1);
    final_filtered(n,2) = filter_values_raw(n,2);
    final_filtered(n,3) = filter_values_raw(n,3); 
end

if final_filtered(n,1) <=0 && final_filtered(n,2) <= 0 &&...
        final_filtered(n,1) <= 0
    final_filtered(n,1) = NaN;
    final_filtered(n,2) = NaN; 
    final_filtered(n,3) = NaN;
end
end

% Eliminate NaNs in data set
final_filtered = final_filtered(0== sum(isnan(final_filtered), 2), :);

% Normalize filtered set
for n=1:length(final_filtered)
sum_final_filtered(n) = final_filtered(n,1) + final_filtered(n,2) + ...
    final_filtered(n,3);
normalized_final_filtered(n,1) = final_filtered(n,1)/sum_final_filtered(n);
normalized_final_filtered(n,2) = final_filtered(n,2)/sum_final_filtered(n);
normalized_final_filtered(n,3) = final_filtered(n,3)/sum_final_filtered(n);
end


%% Compile all region intensity data
    for n=1:length(final_filtered)
   red_f = normalized_final_filtered(n,1);
   green_f = normalized_final_filtered(n,2);
   blue_f = normalized_final_filtered(n,3);
   %Original
    % red-yellow hue
    if red_f >= green_f && green_f >= blue_f
        
        hue_background_final(n) = 60*(green_f-blue_f)/(red_f-blue_f);   
        
    % yellow-green hue
    elseif green_f > red_f && red_f >= blue_f
    
        hue_background_final(n) = 60*(2-(red_f-blue_f)/(green_f-blue_f));    
        
    % Green-Cyan hue    
    elseif green_f >= blue_f && blue_f > red_f
        
       hue_background_final(n) = 60*(2+(blue_f-red_f)/(green_f-red_f));   
        
    % Cyan-Blue hue        
    elseif blue_f > green_f && green_f > red_f
        
     hue_background_final(n) = 60*(4-(green_f-red_f)/(blue_f-red_f));      
    
        % Blue-Magenta hue    
    elseif blue_f > red_f && red_f >= green_f
        
        hue_background_final(n) = 60*(4+(red_f-green_f)/(blue_f-green_f));
           
    % Magenta-Red hue    
    elseif red_f >= blue_f && blue_f > green_f
        
      hue_background_final(n) = 60*(6-(blue_f-green_f)/(red_f-green_f));  
    end
    % 
        % Calculate hues for smoothed background
   red_f15 = rgb_normalized_b15(n,1);
   green_f15 = rgb_normalized_b15(n,2);
   blue_f15 = rgb_normalized_b15(n,3);
   
   red_f30 = rgb_normalized_b30(n,1);
   green_f30 = rgb_normalized_b30(n,2);
   blue_f30 = rgb_normalized_b30(n,3);
      % For filter 15
    % red-yellow hue
    if red_f15 >= green_f15 && green_f15 >= blue_f15
        
        hue_15(n) = 60*(green_f15-blue_f15)/(red_f15-blue_f15);   
        
    % yellow-green hue
    elseif green_f15 > red_f15 && red_f15 >= blue_f15
    
        hue_15(n) = 60*(2-(red_f15-blue_f15)/(green_f15-blue_f15));    
        
    % Green-Cyan hue    
    elseif green_f15 >= blue_f15 && blue_f15 > red_f15
        
       hue_15(n) = 60*(2+(blue_f15-red_f15)/(green_f15-red_f15));   
        
    % Cyan-Blue hue        
    elseif blue_f15 > green_f15 && green_f15 > red_f15
        
     hue_15(n) = 60*(4-(green_f15-red_f15)/(blue_f15-red_f15));      
    
        % Blue-Magenta hue    
    elseif blue_f15 > red_f15 && red_f15 >= green_f15
        
        hue_15(n) = 60*(4+(red_f15-green_f15)/(blue_f15-green_f15));
           
    % Magenta-Red hue    
    elseif red_f15 >= blue_f15 && blue_f15 > green_f15
        
      hue_15(n) = 60*(6-(blue_f15-green_f15)/(red_f15-green_f15));  
    end    
    
   % For filter 30
    % red-yellow hue
    if red_f30 >= green_f30 && green_f30 >= blue_f30
        
        hue_30(n) = 60*(green_f30-blue_f30)/(red_f30-blue_f30);   
        
    % yellow-green hue
    elseif green_f30 > red_f30 && red_f30 >= blue_f30
    
        hue_30(n) = 60*(2-(red_f30-blue_f30)/(green_f30-blue_f30));    
        
    % Green-Cyan hue    
    elseif green_f30 >= blue_f30 && blue_f30 > red_f30
        
       hue_30(n) = 60*(2+(blue_f30-red_f30)/(green_f30-red_f30));   
        
    % Cyan-Blue hue        
    elseif blue_f30 > green_f30 && green_f30 > red_f30
        
     hue_30(n) = 60*(4-(green_f30-red_f30)/(blue_f30-red_f30));      
    
        % Blue-Magenta hue    
    elseif blue_f30 > red_f30 && red_f30 >= green_f30
        
        hue_30(n) = 60*(4+(red_f30-green_f30)/(blue_f30-green_f30));
           
    % Magenta-Red hue    
    elseif red_f30 >= blue_f30 && blue_f30 > green_f30
        
      hue_30(n) = 60*(6-(blue_f30-green_f30)/(red_f30-green_f30));  
    end     
    
    end
    
%% Create plots for analyzing data

figure; scatter3(rgb_normalized(:,1), rgb_normalized(:,2),...
   rgb_normalized(:,3));

figure; scatter3(rgb_raw(:,1), rgb_raw(:,2), rgb_raw(:,3));

figure; hist(Hue, 120)

%%
figure; scatter3(rgb_normalized_filtered(:,1), rgb_normalized_filtered(:,2),...
   rgb_normalized(:,3));

figure; scatter3(rgb_raw_filtered(:,1), rgb_raw_filtered(:,2),...
    rgb_raw_filtered(:,3));

figure; hist(Hue_filtered, 120)

%%
figure; hist(Background_Hue,120)

%% 
figure; hist(hue_background_final,120)
figure; scatter3(normalized_final_filtered(:,1),normalized_final_filtered(:,2),...
    normalized_final_filtered(:,3));

figure; scatter3(final_filtered(:,1),final_filtered(:,2),...
    final_filtered(:,3));

%%
figure; hist(hue_15,360);
figure; hist(hue_30,360);

figure; scatter3(rgb_normalized_b15(:,1),rgb_normalized_b15(:,2),...
    rgb_normalized_b15(:,3));

%% Calculate numbers of particles in each color range
  
% Red

Red_count = sum(hue_background_final > 355 | hue_background_final < 20);
Fuschia_count = sum(hue_background_final > 300 & hue_background_final < 355);
Blue_count = sum(hue_background_final > 231 & hue_background_final < 270);
Light_blue_count = sum(hue_background_final > 181 & hue_background_final < 230);
Blue_green_count = sum(hue_background_final > 135 & hue_background_final < 180);
Yellow_green_count = sum(hue_background_final > 25 & hue_background_final < 110);
Green_count = sum(hue_background_final > 111 & hue_background_final < 130);

sum_count = Red_count+Fuschia_count+Blue_count+Light_blue_count+Blue_green_count...
    +Yellow_green_count+Green_count;

Red_percentage = Red_count/sum_count
Fuschia_percentage = Fuschia_count/sum_count
Blue_percentage = Blue_count/sum_count
Light_blue_percentage = Light_blue_count/sum_count
Blue_green_percentage = Blue_green_count/sum_count
Yellow_green_percentage = Yellow_green_count/sum_count
Green_percentage = Green_count/sum_count

% Calculate intensities
for n=1:length(filter_values_raw)
intensity(n) = (filter_values_raw(n,1) + filter_values_raw(n,2) + ...
    filter_values_raw(n,3))/3;
end

%% Eliminate sections not meeting minimum intensity and/or empty sets
min_filtering_intensity = 0;

for n = 1:length(rgb_raw)
if intensity(n) <= min_filtering_intensity
    final_filtered(n,1) = NaN;
    final_filtered(n,2) = NaN; 
    final_filtered(n,3) = NaN;
else
    final_filtered(n,1) = filter_values_raw(n,1);
    final_filtered(n,2) = filter_values_raw(n,2);
    final_filtered(n,3) = filter_values_raw(n,3); 
end

if final_filtered(n,1) <=0 && final_filtered(n,2) <= 0 &&...
        final_filtered(n,1) <= 0
    final_filtered(n,1) = NaN;
    final_filtered(n,2) = NaN; 
    final_filtered(n,3) = NaN;
end
end

% Eliminate NaNs in data set
final_filtered = final_filtered(0== sum(isnan(final_filtered), 2), :);

% Normalize filtered set
for n=1:length(final_filtered)
sum_final_filtered(n) = final_filtered(n,1) + final_filtered(n,2) + ...
    final_filtered(n,3);
normalized_final_filtered(n,1) = final_filtered(n,1)/sum_final_filtered(n);
normalized_final_filtered(n,2) = final_filtered(n,2)/sum_final_filtered(n);
normalized_final_filtered(n,3) = final_filtered(n,3)/sum_final_filtered(n);
end

%% Filtered histogram


