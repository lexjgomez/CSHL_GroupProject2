%% Directory path across computers
clear; clc;

% list of known computers: point to git repo directory on your device
if exist(fullfile('C:', 'Users', 'kedea', 'Documents', 'CSHL_GroupProject2'),'dir')
    home_dir = 'C:\Users\kedea\Documents\CSHL_GroupProject2';
elseif exist(fullfile('C:', 'Users', 'course', 'Documents', 'CSHL_GroupProject2'),'dir')
    home_dir = 'C:\Users\course\Documents\CSHL_GroupProject2';
elseif exist(fullfile('/', 'Users', 'carolinejia', 'Documents', 'GitHub', 'CSHL_GroupProject2'), 'dir')
    home_dir = fullfile('/Users', 'carolinejia', 'Documents', 'CSHL_GroupProject2');
elseif exist(fullfile('C:', 'Users', 'Admin', 'Documents', 'CSHL_GroupProject2'), 'dir')
    home_dir = fullfile('C:\Users\Admin\Documents\CSHL_GroupProject2');
else
    error('Add your directory to this list (or rewrite this if there is a better way) - Kat')
% code_dir = fullfile('C:', 'Neuda2023', 'Code', 'w1d1', 'npy-matlab-master');
end

addpath(genpath(home_dir))
cd(home_dir)

% check that the data is one this computer in the same exact location (in
% this file to be safe - git will ignore data (.gitignore) so it needs to
% be added to the repo on your device
if ~exist('allData','dir')
    error('Please add the allData folder from Steinmetz into this repo - Kat')
end

%% Load data
whichMouse = 'Cori_2016-12-14';

% genpath now sees all of the data inside allData
sesPath = fullfile(home_dir, 'allData', whichMouse);
[S, neurons, trials, regions] = openSession(sesPath);  

%% Select one area in animal
region = 'VISp';
% get region index
region_index = find(strcmp(regions.name, region));
% identify neurons in region
region_neurons = find(neurons.region == region_index);

% How many neurons in each region?
num_unit_all = length(region_neurons);
num_sample = 51;
num_trial = trials.N;

unit_selected = randsample(length(region_neurons), num_sample);

%% Clean up data
%get rid of outliers

% pick out units that responded after stimulus, 99% confidence


%% get spike data and downsample
spike_counts = get_spike_counts(trials, S, region_neurons);
spike_counts_downsample = spike_counts(:, unit_selected);


%% Multiple regression

Predicted = imultipleregress(spike_counts_downsample);

% R-squared calculation
RSS = (spike_counts_downsample - Predicted) .^ 2;
RSS = mean(RSS, 1);
TSS = (spike_counts_downsample - mean(spike_counts_downsample, 1)) .^ 2;
TSS = mean(TSS, 1);

R_squared = 1 - (RSS./TSS);
figure; histogram(R_squared);
xlabel = 'R'


%% Fano Factors 

% example cell figure 

% choose high predictor cells 

%% Multiple Regression Behavior from Neurons



%% Scaling up - apply to multiple brain areas 
