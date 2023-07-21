function [dat, data] = loadsession(path2data, fishID)
%% Open and QC and display one zebrafish data set 
% Reads data downloaded from https://figshare.com/articles/dataset/Whole-brain_light-sheet_imaging_data/7272617
% First set fishID: fish 4 is good to start with. 
if ~exist('fishID', 'var'), fishID = 'subject_4'; end
% fpsec = frames per second

%% load accessories in .mat file - this section and next take ~10 sec
infoFile = fullfile(path2data, fishID, 'data_full.mat');
load(infoFile) % load information
times = (1:length(data.stim_full))'/data.fpsec; % vector of imaging times

validCells = setdiff( 1:size(data.CellXYZ,1), data.IX_inval_anat);
N = length(validCells);
data.CellXYZ = data.CellXYZ(validCells,:);
data.CellXYZ_norm = data.CellXYZ_norm(validCells,:);
%plot3(data.CellXYZ_norm(:,1),data.CellXYZ_norm(:,2),data.CellXYZ_norm(:,3),'.') % z-values inferred from imaging planes
%xlabel('X'); ylabel('Y'); zlabel('Z'); title('3D Locations of Cells'); grid


%% load calcium traces in HDF5 file and trim suspect frames - 5-10 sec - 30 sec if display 
datFile = fullfile(path2data, fishID, 'TimeSeries.h5');% data file
dat = h5read(datFile,'/CellResp'); % ,[1 1],[15000 1650] for subset

% The last two frames of fish 7 time series seem discontinuous with the preceding images
if fishID == "subject_7"
    goodFrames = 1:size(dat,2)-2;
else
    goodFrames = 1:size(dat,2);
end
if length(goodFrames) < size(dat,2) % 
    times=times(goodFrames);
    dat = dat(goodFrames,:);
    data.Behavior_full = data.Behavior_full(:,goodFrames);
    data.stim_full = data.stim_full(goodFrames);
    data.Eye_full = data.Eye_full(:,goodFrames);
end
numFrames = size(dat,1);
%
%figure; imagesc(times,1:size(dat,1), sqrt(dat)); grid % most values around 1; skewed distn with values up to ~14
fluomeans = mean(dat,2);
fluoSDs = std(dat,0,2);


%% If not already log-scale, then Log transform and transpose
% I'm checking crudely here...
% Some transform should be done for data from each fish depending on characteristics
if all(dat(:) > 0)
    dat = log(dat'); "Log transform"
else
    dat = dat';
end

%% Crude drift correction
% should check for drift in a more sophisticated way
xx = [ ones(size(dat,1),1) times]; % predictors for drift
tmp = xx \ dat; % estimate mean drift rates
drifts = tmp(2,:)'; intcpts = tmp(1,:)';
% typically offsets between 0 and 0.5 for fish 1
% typical drifts between -5E-4 and 5E-4 for fish 1
% these are actually substantial: a difference of 0.5 for 1250 cells and
% greater than half of cell SD for 24K neurons

dat = dat - xx * tmp ; % compensate drift



