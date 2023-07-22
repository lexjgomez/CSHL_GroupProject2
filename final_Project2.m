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
%num_sample = 51;

rng(10); %setting our random number generator
%unit_selected = randsample(length(region_neurons), num_sample); %randomly selecting neurons

%% get spike data 

[spike_counts, allSpikes, allSpikesperTrial] = get_spike_counts(trials, S, region_neurons); %neurons x times (in 100 ms bins, 5 pre, 6 post trial onset) x num trial


%downsample everything
% spike_counts_downsample = spike_counts(:, unit_selected);
% allSpikes_downsample = allSpikes(unit_selected,:);
% allSpikesperTrial_downsample = allSpikesperTrial(unit_selected,:,:);

%% Plot some raw data 

figure;
imagesc(sum(allSpikesperTrial_downsample,3))
colorbar
xlabel('Time Bins [100 ms]')
ylabel('Neuron')
title('Summed Spikes Over Time Bins per Neuron')


%% Plot neuron by neuron 



%% Pull Out Responsive Cells and Downsample

% pick out units that responded after stimulus, paired t-test
num_unit_all = length(region_neurons);
[responsiveUnits, nonresponsiveUnits, Indices] = classify_by_response(allSpikesperTrial,num_unit_all);

%extract indices
respInd = Indices.Trial_IndiciesResponsive;
nonrespInd = Indices.Trial_IndiciesNonResponsive;

%get spike counts at these indices
spike_counts_resp = spike_counts(:,respInd);
spike_counts_nonresp =spike_counts(:,nonrespInd);

num_sample = 51;

% select random subsample of units 
unit_selected_resp = randsample(size(responsiveUnits,1), num_sample);
unit_selected_nonresp = randsample(size(nonresponsiveUnits,1), num_sample);

%get responsive and non-responsive by neuron x bins x trials
selected_Responsive = responsiveUnits(unit_selected_resp,:,:);
selected_Unresponsive = nonresponsiveUnits(unit_selected_nonresp,:,:);

%get responsive and non-responsive by times x neurons THIS IS WHAT YOU USE
%FOR THE REGRESSION ANALYSIS
spike_counts_resp_downsample = spike_counts_resp(:,unit_selected_resp);
spike_counts_nonresp_downsample = spike_counts_nonresp(:,unit_selected_nonresp);


%normalize for plotting 
norm_Selected_Responsive = selected_Responsive./max(selected_Responsive,2);
norm_Selected_Unresponsive = selected_Unresponsive./max(selected_Unresponsive,2);


s = figure
s.Position = [776.2,524.2,1031.2,420];
subplot(1,2,1)
plot(squeeze(mean(norm_Selected_Responsive,1)),'r')
hold on
plot(mean(squeeze(mean(norm_Selected_Responsive,1)),2),'k', 'LineWidth',3)
xlabel('bins')
ylabel('normalized firing rate')
title('Subsample of Neurons Classified as *Responsive*')

subplot(1,2,2)
plot(squeeze(mean(norm_Selected_Unresponsive,1)),'b')
hold on
plot(mean(squeeze(mean(norm_Selected_Unresponsive,1)),2),'k', 'LineWidth',3)
xlabel('bins')
ylabel('normalized firing rate')
title('Subsample of Neurons Classified as *NON-Responsive*')

% spike_counts_downsample = responsiveUnits(unit_selected_resp,:,:);


%% Multiple regression

%for responsive
Predicted = imultipleregress(spike_counts_resp_downsample);

% R-squared calculation
RSS = (spike_counts_resp_downsample - Predicted) .^ 2;
RSS = mean(RSS, 1);
TSS = (spike_counts_resp_downsample - mean(spike_counts_resp_downsample, 1)) .^ 2;
TSS = mean(TSS, 1);

R_squared = 1 - (RSS./TSS);

%for unresponsive
Predicted_non = imultipleregress(spike_counts_nonresp_downsample);

% R-squared calculation
RSS_n = (spike_counts_nonresp_downsample - Predicted_non) .^ 2;
RSS_n = mean(RSS_n, 1);
TSS_n = (spike_counts_nonresp_downsample - mean(spike_counts_nonresp_downsample, 1)) .^ 2;
TSS_n = mean(TSS_n, 1);

R_squared_n = 1 - (RSS_n./TSS_n);
figure 
subplot(1,2,2)
histogram(R_squared_n);
xlabel ( 'Non-responsive R²')

%% Fano Factors 

%spike_counts_downsample = time x neurons 
% example cell figure 
fanos = ifanofactor (allSpikes_downsample);
% choose high predictor cells 

%% Multiple Regression Behavior from Neurons



%% Scaling up - apply to multiple brain areas 
