%%% Script for determining HUVEC coverage in planar vasc networks
%%% Ian Kinstlinger & Kristen Means
%%% MillerLab
%%% rev. 11/2021
%%% Cell shape for circularity


% Dependency: must have uipickfiles.m on your active path
% Interactive ROI functions require MATLAB 2019b 

imList = uipickfiles('out','struct','Prompt','Select images');
%Select all the images to process (not the directory but the image files)

for i = 1:length(imList)

im = im2double(imread(imList(i).name));

if mean(mean(im)) < .02 
%     im = 2^4.*im; %fix 12 bit images which Matlab reads as 16bit
    %comment this out or adjust accordingly if the souce images are not 12-bit!
    
end

% Adjust data to span data range.
im = imadjust(im);

% Threshold image - adaptive threshold (sensitivity normally 0.7)
bin = imbinarize(im, 'adaptive', 'Sensitivity', 0.7000, 'ForegroundPolarity', 'bright');

%%Define an ROI large enough to enclose the in focus monolayer
roiX = 800; roiY = 800; 
roiBox = figure(); imshow(im);

% Draw an ROI of the size specified on each image, but let the user move
% the ROI to cover the desired region
roi = drawrectangle('Position', [1024-(roiX/2), 1024-(roiY/2),roiX,roiY], 'InteractionsAllowed', 'translate');

%Double-click final ROI to submit
wait(roi) 
roipos = roi.Position;
close(roiBox)

%Crop original image to ROI area
cropOriginal = imcrop(im, roipos); 
figure(); imshow(cropOriginal);

%Crop binary image to ROI area
crop = imcrop(bin, roipos); 
figure(); imshow(crop);

%Clean noise smaller than XX pixels
cropCleaned = bwareaopen(crop,500);
%imshow(cropCleaned)

% Use distance transform to create catchment basins 
D = -bwdist(~cropCleaned);
%figure(); imshow(D,[])

% Perform watershed segmentation on catchment basins
Ld = watershed(D);
%figure(); imshow(label2rgb(Ld))

% Filter out tiny local minima 
mask = imextendedmin(D,2);
%figure(); imshowpair(cropCleaned,mask,'blend')

% Modify the distance transform so that no minima occur at the filtered-out locations
D2 = imimposemin(D,mask);
Ld2 = watershed(D2);
cropCleanedSeg = cropCleaned;
cropCleanedSeg(Ld2 == 0) = 0;
%figure(); imshow(cropCleanedSeg)

%Fill holes in segmented image
cropCleanedSegFill = imfill(cropCleanedSeg,'holes');
%figure(); imshow(cropCleanedSegFill)

% Create masked image of cropped ROI.
maskedImageCrop = cropOriginal;
maskedImageCrop(~cropCleanedSegFill) = 0;
figure(); imshow(maskedImageCrop);

% Display segmented cells with colormap and boundaries
[B,L] = bwboundaries(cropCleanedSegFill,'noholes');
figure(); imshow(label2rgb(L,@jet,[.5 .5 .5]))
hold on
for k = 1:length(B)
  boundary = B{k};
  plot(boundary(:,2),boundary(:,1),'w','LineWidth',2)
end

% Make structure of cell properties
stats = regionprops(L,'Area','Centroid','Perimeter','Solidity','Circularity','MajorAxisLength','MinorAxisLength','MaxFeretProperties','MinFeretProperties');

% Loop over the boundaries
for k = 1:length(B)

  % Obtain (X,Y) boundary coordinates corresponding to label 'k'
  boundary = B{k};

  % Obtain the perimeter & solidity calculation corresponding to label 'k'
  perimeter = stats(k).Perimeter;
  solidity = stats(k).Solidity;
  
  % Obtain the area calculation corresponding to label 'k'
  area = stats(k).Area;
  
  % Compute the roundness metric
  metric = 4*pi*area/perimeter^2;
  
  % Display the results
  metric_string = sprintf('%2.2f',area);
  text(boundary(1,2)-35,boundary(1,1)+13,metric_string,'Color','black',...
       'FontSize',10,'FontWeight','bold')
  
end


%----------
% New cropped area with different threshold sensitivity for cell coverage
% calculation (typically higher sensititivity to catch more of the signal, 
% note individual cells are not segmented as well)
% All other watershed steps are the same as above, but cell coverage area is
% calculated by the fraction of GFP pixels over all pixels in selected ROI

% Threshold image - adaptive threshold (sensitivity normally 0.9)
bin = imbinarize(im, 'adaptive', 'Sensitivity', 0.9000, 'ForegroundPolarity', 'bright');

%%Define an ROI large enough to enclose the in focus monolayer
roiX = 1800; roiY = 600; 
roiBox = figure(); imshow(im);

% Draw an ROI of the size specified on each image, but let the user move
% the ROI to cover the desired region
roi = drawrectangle('Position', [1024-(roiX/2), 1024-(roiY/2),roiX,roiY], 'InteractionsAllowed', 'translate');
wait(roi) %Double-click final ROI to submit
roipos = roi.Position;
close(roiBox)

%Crop original image to ROI area
cropOriginal = imcrop(im, roipos); 
figure(); imshow(cropOriginal);

%Crop binary image to ROI area
crop = imcrop(bin, roipos); 
figure(); imshow(crop);

%Clean noise smaller than XX pixels
cropCleaned = bwareaopen(crop,500);
%imshow(cropCleaned)

% Use distance transform to create catchment basins 
D = -bwdist(~cropCleaned);
%figure(); imshow(D,[])

% Perform watershed segmentation on catchment basins
Ld = watershed(D);
%figure(); imshow(label2rgb(Ld))

% Filter out tiny local minima 
mask = imextendedmin(D,2);
%figure(); imshowpair(cropCleaned,mask,'blend')

% Modify the distance transform so that no minima occur at the filtered-out locations
D2 = imimposemin(D,mask);
Ld2 = watershed(D2);
cropCleanedSeg = cropCleaned;
cropCleanedSeg(Ld2 == 0) = 0;
%figure(); imshow(cropCleanedSeg)

%Fill holes in segmented image
cropCleanedSegFill = imfill(cropCleanedSeg,'holes');
%figure(); imshow(cropCleanedSegFill)  

% Create masked image of cropped ROI.
maskedImageCrop = cropOriginal;
maskedImageCrop(~cropCleanedSegFill) = 0;
figure(); imshow(maskedImageCrop);

% Find cell coverage from segmented image 
totalArea = roiX*roiY;
sumPixels = sum(sum(cropCleanedSegFill)); %Count all signal-positive pixels
frac = sumPixels/totalArea %Report result as a fraction of total ROI area

%Add coverage field to structure stats 
stats(1).Coverage = frac;

% % *Export data* - uncomment and name file when ready to save 
%   File will save to the current folder as an excel sheet, Warning: files
%   will be saved over if the name is not changed between runs

%filename = 'testFilename.xlsx';
%writetable(struct2table(stats),filename,'Sheet',1,'Range','D1')

end

