

function [trackStats, spots] = postprocessSimpleTracks(inputFile, ...
    frameGap, numFrames, saveExcelFile, makeGraphs)

% POSTPROCESSSIMPLETRAKCS Computes tracking stats from out .csv file
% produced by TrackMate.
%
% INPUT inputFile: The located on the .csv file produced by TrackMate                     
%       frameGap: time between frames (in seconds)                   
%       numFrames: total number of frames in movie
%       saveExcelFile: if true will save the output table as an excel file              
%       makeGraphs: if true will produce graphs of spot count of time
%
% OUTPUT  trackStats: struct containing track statistics
%         spots: struct containing spot statistics
% 
% REMARKS: This function processes the outrput file from TrackMate produced
% by selecting "Export all spot statistics". It includes tracks with only a
% single spot and does not account for gap linking or track splitting and
% merging. For more information on tracking see: 
% Tinevez, Jean-Yves, et al. Methods 115 (2017): 80-90.
%
%
% Author: Jeremy Pike
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Load excel file
spreadsheetNum = xlsread(inputFile);

% Pull out spot track IDs, frames, quality and mean intensity
spots.trackID = spreadsheetNum(:, 4) + 1;
spots.frame = spreadsheetNum(:, 10) + 1;
spots.quality = spreadsheetNum(:, 5);
spots.meanIntensity = spreadsheetNum(:, 14);

% Find total number of spots and number of tracks
numSpots = size(spreadsheetNum, 1);
numTracks = max(spots.trackID);

% Identify location of isolated spots not in tracks
isolatedSpots = isnan(spots.trackID);


% Assign all isolated spots a new track number
numTracksInc = numTracks;
for i = 1: numSpots
    if isolatedSpots(i) == true
         numTracksInc = numTracksInc + 1;
        spots.trackID(i) = numTracksInc;    
    end
end

% Create emtpy arrays to gtoup spots into tracks
tracks = [];
for t = 1 : numTracksInc
    tracks{t, 1}.spotID = [];
    tracks{t, 1}.frame = [];
    tracks{t, 1}.quality = [];
    tracks{t, 1}.meanIntensity = [];
end

% Group spots and their properties into tracks
for i = 1: numSpots
    tracks{spots.trackID(i, 1)}.spotID = ...
        [tracks{spots.trackID(i, 1)}.spotID; i];
    tracks{spots.trackID(i, 1)}.frame = ... 
        [tracks{spots.trackID(i, 1)}.frame; spots.frame(i)];
    tracks{spots.trackID(i, 1)}.quality = ... 
        [tracks{spots.trackID(i, 1)}.quality; spots.quality(i)];
    tracks{spots.trackID(i, 1)}.meanIntensity = ... 
       [tracks{spots.trackID(i, 1)}.meanIntensity; spots.meanIntensity(i)];
 
end

% Calculate key statistics for each track
trackStats = [];
for t = 1 : numTracksInc
    trackStats.numSpots(t, 1) = size(tracks{t,1}.spotID, 1);
    trackStats.duration(t, 1) ... 
        =  (trackStats.numSpots(t, 1) - 1) * frameGap;
    trackStats.meanQuality(t, 1) = mean(tracks{t,1}.quality, 1);
    trackStats.meanMeanIntensity(t, 1) ... 
        = mean(tracks{t,1}.meanIntensity, 1);
    trackStats.firstFrame(t,1) = min(tracks{t,1}.frame, [], 1);
    trackStats.lastFrame(t,1) = max(tracks{t,1}.frame, [], 1);
end

% Make a table from trackStats
trackStatsTable = table(trackStats.numSpots, trackStats.duration, ... 
    trackStats.meanQuality, trackStats.meanMeanIntensity, ... 
    trackStats.firstFrame, trackStats.lastFrame);
trackStatsTable.Properties.VariableNames = {'numSpots' 'duration' ...
    'meanQuality' 'meanMeanIntensity' 'firstFrame' 'lastFrame'};  

% If specified save table as an excel file
if saveExcelFile
    writetable(trackStatsTable, [inputFile(1 : length(inputFile) - 4) ...
        '_trackStats.xls']) 
end

% Calculate number of spots per frame
spotsPerFrame = zeros(numFrames, 1);
spotPerFrameCum = zeros(numFrames, 1);
for i = 1 : numSpots
    spotsPerFrame(spots.frame(i), 1) = ... 
        spotsPerFrame(spots.frame(i), 1) + 1;
end
% Calculate cumulative spots per frame 
for t = 1 : numFrames
    spotPerFrameCum(t, 1) = sum(spotsPerFrame(1 : t, 1));
end

% If specified plot graphs
if makeGraphs
    figure; plot(1:1:numFrames, spotsPerFrame)
    xlabel('frame number')
    ylabel('spots per frame')
    figure; plot(1:1:numFrames, spotPerFrameCum)
    xlabel('frame number')
    ylabel('spots per frame (cumulative)')
end