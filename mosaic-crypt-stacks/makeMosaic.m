
function makeMosaic(minSizeDiffPercent, minDistance, saveMosaics)

% MAKEMOSAIC Creates a mosaic figure of clone aquisitions from search and 
% find confocal scan. All .lif files in a user specifed directory are 
% processed. Duplicates clone aquistions are removed and displayed in a 
% seperate mosaic figure. 
%
% INPUT minSizeDiffPercent: minimum percentage area difference between
%                           clones to be classed as not duplicates
%       minDistance: minimum distance (microns) between clones to be
%                    classed as not duplicates
%       saveMosaics: boolean, if true will save the mosaics to specified
%                    directory
%
% REMARKS: This is a custom function for use by Doug Wintons lab at the
% CRUK-CI
%
% DEPENDENCIES: Requires the Bioformats MATLAB toolbox which can be
% downloaded from http://downloads.openmicroscopy.org/bio-formats. Version
% 5.3.4 was used to test this script.
%
% Author: Jeremy Pike
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set default input parameters if not provided
if nargin < 1
    minSizeDiffPercent = 10;
end
if nargin < 2
    minDistance = 50;
end
if nargin < 3
    saveMosaics = true;
end

%% Load data and perform maximal projection
display('Loading data ....')
% Ask user to specify a directory containing data
folderPath = uigetdir('','Specify input data folder');
% Ask user to specify a directory for saving the output tiled images
outputFolderPath = uigetdir('','Specify output folder');
% Set current folder to specified directory
cd(folderPath);
% Find all leica .lif files in the specified directory
leicaFiles = dir('*.lif');

if (isempty(leicaFiles))
   display('Please provide a directory with .lif files')
   return
end

nFiles = length(leicaFiles);


for i = 1 : nFiles
   % retrieve path for current file
   currentfilename = leicaFiles(i).name;
   
   % load and poject all series in the file, also retrieve stage
   % coordinates
   [seriesDataMaxProj, seriesPlanePosition]  ... 
       = bioformReadMaxProj(currentfilename);
   
   % add file data to cell arrays for all files
   if i > 1
       dataMaxProj = [dataMaxProj; seriesDataMaxProj];
       planePositions = [planePositions; seriesPlanePosition];
   else 
       dataMaxProj = seriesDataMaxProj;
       planePositions = seriesPlanePosition;
   end
   
   display(['File ' num2str(i) ' of ' num2str(nFiles) ' imported'])
end

% total number of series across all files
imageCount = size(dataMaxProj, 1);


%% Calculate sizes
display('Calculating clone sizes and sorting ....')

% create vector counting all data from the second channel across all series
dataAllVec = [];
for i = 1 : imageCount
    tempVec = dataMaxProj{i, 1}(:, :, 2);
    tempVec = tempVec(:);
    dataAllVec = [dataAllVec; tempVec];
end
% calculate otsu threshold across all data series
otsuLevel = multithresh(dataAllVec);
% clear the all data vector from memory as no longer needed
clear dataAllVec tempVec

% calculate clone sizes using threshold
cloneSizes = zeros(imageCount, 1);
for i = 1 : imageCount
        
        % Apply threshold to generate binary image
        BW = dataMaxProj{i, 1}(:, :, 2) >= otsuLevel;  
        % Calculate areas of all connected components in binary iamge
        R = regionprops(bwconncomp(BW), 'Area');
        % Find the largest area
        if isempty(cat(1, R.Area)) == false
            
            cloneSizes(i, 1) = max(cat(1, R.Area));
        end
            
end

% find indexes for sorting by size in descending order
[~, orderedIndex] = sort(cloneSizes, 'descend');

%% Find clone duplicates
display('Finding duplicates ....')
duplicate = zeros(imageCount, 1);

% first clone is never a duplicate
duplicate(i, 1) = false;
% Calculate distances and size difference between clones
for i = 2 : imageCount
    % only need to calculate to the ith clone so that first instance of
    % clone is not registered as a duplcicate
    sizeDifPer = zeros(i - 1, 1);
    distances = zeros(i - 1, 1);
    for j = 1 : i - 1
    
       % percent size difference
       sizeDifPer(j, 1) = abs(cloneSizes(i) - cloneSizes(j)) ... 
           / cloneSizes(i);
       % distance between clones based on stage coordinates
       distances(j, 1) = sqrt(( ... 
           planePositions(i, 1) - planePositions(j, 1))^2 + ... 
           (planePositions(i, 2) - planePositions(j, 2))^2);
 
   end
    % true if minimum percentage size difference less than threshold
    sizeCheck = min(sizeDifPer) < minSizeDiffPercent / 100;
    % true if minimum distance less than threshold
    distanceCheck = min(distances) < minDistance;
    % if both checks are true then clone is a duplicate
    if (sizeCheck && distanceCheck)
        duplicate(i, 1) = true;
    else
        duplicate(i, 1) = false;
    end
   
end

    
%% Create mosaics
display('Making mosaics ...')

% for non-duplicate clones
figure
mosaicNonDup = ...
    mosaicFigure(dataMaxProj, duplicate == false, orderedIndex, 20);
% for duplicates
figure
mosaicDup = mosaicFigure(dataMaxProj, duplicate, orderedIndex, 20);

% find position of final '\' in input directory for output file names
slashes = strfind(folderPath, '\');
% if requsted save mosaic images 
if saveMosaics
    imwrite(mosaicNonDup,[outputFolderPath '\mosaicData_' ... 
        folderPath(slashes(length(slashes)) + 1 : end) '.tif'], 'tif') 
    imwrite(mosaicDup,[outputFolderPath '\mosaicDuplicates_' ... 
        folderPath(slashes(length(slashes)) + 1 : end) '.tif'], 'tif') 
end

display('All done :)')
display(['There where ' num2str(sum(duplicate == false)) ... 
    ' non-duplicates and ' num2str(sum(duplicate)) ' duplicates.'])

end

function mosaic = mosaicFigure(data, displayCheck, sortingIndexes, ...
    gapWidth)

% MOASAICFIGURE creates the mosaic figure in a specified order
%
% INPUT data: cell array containing the data from all series
%       displayCheck: boolean vector, where a true value indicates that 
%                     the series should be included in the figure
%       sortingIndexes: vector defining the order for each series in the 
%                       mosaic
%       gapWidth: gap between images in mosaic (pixels)
%
% OUTPUT moasic: 2D array containing the mosaic image
%
% REMARKS: 
%
% Author: Jeremy Pike
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the dimension of each series image
tileRows = size(data{1, 1}, 1);
tileCols =size(data{1, 1}, 2);    

% find the total number of images to be shown
numDisplayImages = sum(displayCheck);

% find number of rows an columns the mosaic image should have 
mosaicCols = ceil(sqrt(numDisplayImages));
mosaicRows = ceil(numDisplayImages / mosaicCols);

% create empty matrix for mosaic
mosaic = zeros(mosaicRows * tileRows + (mosaicRows - 1) * gapWidth, ...
mosaicCols * tileCols + (mosaicCols - 1) * gapWidth, 3);

% keeps track of how many images have been added to mosaic
displayCount = 0;
% loop through all series
for i=1:size(data,1);
   
   % if the current series (specified by sorting index) is to be shown 
   if displayCheck(sortingIndexes(i),1) == true
       
       displayCount = displayCount + 1;
       % find row and column position
        row = ceil(displayCount / mosaicCols);
        col = rem(displayCount, mosaicCols);
        if col == 0
            col = mosaicCols;
        end
        % add data to mosaic in correct position
        mosaic(1 + (row - 1) * (tileRows + gapWidth) ... 
            : tileRows + (row - 1) * (tileRows + gapWidth), ... 
            1 + (col - 1) * (tileCols + gapWidth) : tileCols +(col - 1) ... 
            * (tileCols + gapWidth), :)= data{sortingIndexes(i), 1};
   end
end

% scale each channel between 0 and 1
for i=1:3
    channel = mosaic(:, :, i);
    mosaic(:,:,i) = (channel - min(channel(:))) ... 
        / (max(channel(:)) - min(channel(:)));
end

% swap the red and blue channels for display purposes
temp=mosaic(:,:,1);
mosaic(:, :, 1)=mosaic(:, :, 2);
mosaic(:, :, 2)=temp;

% show the mosaic
imshow(mosaic)


end

