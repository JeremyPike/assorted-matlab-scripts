function [dataSeriesMaxProj, planePosition]  = bioformReadMaxProj(filePath)

% BIOFORMREADMAXPROJ using a Bioformats reader to load all series in a
% specified file and then return the maximal projection of each series. The
% stage coordinates of each series are also returned.
%
% INPUT  filePath: path for file containing image data to be read by
%                  bioformats
%
% OUTPUT dataSeriesMaxProj: cell array where each element contains the
%                           maximally prjected data of each series
%        planePosition: array containing the stage coordinates for each 
%                       series (in microns)
%
% DEPENDENCIES: Requires the Bioformats MATLAB toolbox which can be
% downloaded from http://downloads.openmicroscopy.org/bio-formats. Version
% 5.3.4 was used to test this script.
%
% REMARKS: 
%
% Author: Jeremy Pike
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get a Bioformats reader for the specified file
reader = bfGetReader(filePath);
% Retrieve the number of series from the reader
numSeries = reader.getSeriesCount();
% Retrieve the file metadata
omeMeta = reader.getMetadataStore();

% Empty array for stage coordinates
planePosition = zeros(numSeries, 2);

for i = 1 : numSeries
    % Set the reader to the current series
    reader.setSeries(i - 1);
    
	% Get data dimensions for current series
    rows = omeMeta.getPixelsSizeY(i - 1).getValue();
    cols = omeMeta.getPixelsSizeX(i - 1).getValue();
    numSlices = omeMeta.getPixelsSizeZ(i - 1).getValue();
    numTimePoints = omeMeta.getPixelsSizeT(i - 1).getValue(); 
    numChannels = omeMeta.getPixelsSizeC(i - 1).getValue(); 
    
    % Get stage coordiantes (convert to microns)
    planePosition(i, 1) = omeMeta.getPlanePositionX(i - 1, 0).value() ... 
                        .doubleValue() * 10^6;
    planePosition(i, 2) = omeMeta.getPlanePositionY(i - 1, 0).value() ... 
                        .doubleValue() * 10^6;
    
    % Empty array for series data
    dataSeries = zeros(rows, cols, numSlices, numChannels, numTimePoints);
    
    % For every image plane
    for c = 1 : numChannels
        for t = 1 : numTimePoints
            for z = 1 : numSlices
                % Get plane index for current channel, timepoint and slice
                iPlane = reader.getIndex(z - 1, c -1, t - 1) + 1;
                % Load the plane
                dataSeries(:, :, z, c, t) = bfGetPlane(reader, iPlane);
            end
        end
    end
    
    % Calcaulate the axial maximal projection 
    dataSeriesMaxProj{i, 1} = squeeze(max(dataSeries, [], 3));
    
end

end

