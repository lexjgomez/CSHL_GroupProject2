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

% set random number generator seed 
rng(5); 

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


%% get spike data 

% neurons x times (in 100 ms bins, 5 pre, 6 post trial onset) x num trial
[spike_counts, allSpikes, allSpikesperTrial] = get_spike_counts(trials, S, region_neurons); 

%% PSTH gui plotting
spike_PSTH = get_spike_PSTH(trials, S, region_neurons);
PSTH_GUI(spike_PSTH, trials);


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


s = figure;
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


figure;
subplot(1,2,1)
imagesc(sum(selected_Responsive,3))
colorbar
xlabel('Time Bins [100 ms]')
ylabel('Neuron')
title('Summed Responsive Spikes')
subplot(1,2,2)
imagesc(sum(selected_Unresponsive,3))
colorbar
xlabel('Time Bins [100 ms]')
ylabel('Neuron')
title('Summed Unresponsive Spikes')


%% Multiple Regression - ALL TIMES

%for unresponsive
Predicted_non = imultipleregress(spike_counts_nonresp_downsample);
% R-squared calculation
RSS_n = (spike_counts_nonresp_downsample - Predicted_non) .^ 2;
RSS_n = mean(RSS_n, 1);
TSS_n = (spike_counts_nonresp_downsample - mean(spike_counts_nonresp_downsample, 1)) .^ 2;
TSS_n = mean(TSS_n, 1);

unr_R_squared = 1 - (RSS_n./TSS_n);


%for responsive
Predicted = imultipleregress(spike_counts_resp_downsample);
% R-squared calculation
RSS = (spike_counts_resp_downsample - Predicted) .^ 2;
RSS = mean(RSS, 1);
TSS = (spike_counts_resp_downsample - mean(spike_counts_resp_downsample, 1)) .^ 2;
TSS = mean(TSS, 1);

resp_R_squared = 1 - (RSS./TSS);

figure 
subplot(1,2,1)
histogram(resp_R_squared, -0.1:0.1:0.6);
xlabel('Responsive R²')
ylim([0 35])
set(gca, 'box', 'off');
set(gca, 'tickdir', 'out');

subplot(1,2,2)
histogram(unr_R_squared, -0.1:0.1:0.6);
xlabel('Non-responsive R²')
ylim([0 35])
set(gca, 'box', 'off');
set(gca, 'tickdir', 'out');
linkaxes

%% Perform Multiple Regression on varying number of neurons 
%this takes forever to run
R_squared = iVarymultipleregress(spike_counts_resp_downsample);

avgR2 = mean(R_squared,1);
err = std(R_squared)/sqrt(length(R_squared));

figure()
errorbar([1:num_sample-1], avgR2, err)
ylim([0,0.4])
xlabel('Number of Predictor Neurons'); ylabel('R²');
title(" The effect of increasing the number of Predictors")


%% Fano Factors 

%spike_counts_downsample = time x neurons 
% figures for responsive and unresponsive cell groups
resp_fanos = ifanofactor (selected_Responsive);
unr_fanos = ifanofactor (selected_Unresponsive);


%% Fano factors compared to r^2 values

% boxplot of responsive and non responsive cells split by r^2 and fano
% factors
figure; boxplot([resp_fanos, resp_R_squared', unr_fanos, unr_R_squared'])
xticklabels({'Responsive FFs', 'Responsive r^2', 'Unrsponsive FFs', 'Unresponsive r^2'})
title('Overview of results by responsiveness')

% correlation value
resp_cor = corr(resp_fanos,resp_R_squared');
cor2 = corr(unr_fanos,unr_R_squared');

% correlation figure 
figure; 
title('Relationship between Fano Factor and R^2')

subplot(1,2,1)
scatter(resp_fanos,resp_R_squared)
legend(['r = ' num2str(resp_cor)]);
ylabel('R squared')
xlabel('Fano Factor')
subtitle('Reponsive Neurons')

subplot(1,2,2)
scatter(unr_fanos,unr_R_squared)
legend(['r = ' num2str(cor2)]);
ylabel('R squared')
xlabel('Fano Factor')
subtitle('Unreponsive Neurons')


%% Multiple Regression Part DEUX - Pre and Post Stimulus

%get the stimulus times in the same format as spike_counts
stim_counts = get_trial_counts(trials);

%concatenate pre times and post times

concatResponsivePre = []; concatResponsivePost = [];
concatUnresponsivePre = []; concatUnresponsivePost = [];
for stimIdx = 1:length(stim_counts)
    curStimCount = stim_counts(stimIdx);
    if curStimCount == 1
        %get time window of stim
        % curPreTimesResp =  spike_counts_resp_downsample((stimIdx-5):(stimIdx+4));
        
        curPreTimesResp = spike_counts_resp_downsample(((stimIdx-5):stimIdx-1),:);
        curPostTimesResp = spike_counts_resp_downsample((stimIdx:(stimIdx+4)),:);
            
            concatResponsivePre = [concatResponsivePre; curPreTimesResp];
            concatResponsivePost = [concatResponsivePost; curPostTimesResp];

        curPreTimesNonResp = spike_counts_resp_downsample(((stimIdx-5):stimIdx-1),:);
        curPostTimesNonResp = spike_counts_resp_downsample((stimIdx:(stimIdx+4)),:);
                
            concatUnresponsivePre = [concatUnresponsivePre; curPreTimesNonResp];
            concatUnresponsivePost = [concatUnresponsivePost; curPostTimesNonResp];

    end 
end 

all_concats = zeros(4,size(concatUnresponsivePre,1),size(concatUnresponsivePre,2));
all_concats(1,:,:) = concatResponsivePre;
all_concats(2,:,:) = concatResponsivePost;
all_concats(3,:,:) = concatUnresponsivePre;
all_concats(4,:,:) = concatUnresponsivePost;
all_concat_names = {'Resp_Pre','Resp_Post','Unresp_Pre','Unresp_Post'};

numTrials = length(trials.visStimTime);

for feedMeSeymour = 1:4

    curConcats = squeeze(all_concats(feedMeSeymour,:,:));
    curConcatName = char(all_concat_names(feedMeSeymour)); 


    %for responsive
    Predicted = imultipleregress(curConcats);
    
    % R-squared calculation
    RSS = (curConcats - Predicted) .^ 2;
    RSS = mean(RSS, 1);
    TSS = (curConcats - mean(curConcats, 1)) .^ 2;
    TSS = mean(TSS, 1);
    
    concat_R_squared = 1 - (RSS./TSS);
    
    figure(200)
    subplot(1,4,feedMeSeymour)
    histogram(concat_R_squared);
    xline(0,'r')
    xlabel([curConcatName ' R²'])
    title(curConcatName)

    %fano factors
    %reshape shit so that the fanny factors work
    concatFannyShaped = [];
    for fannyShape = 1:num_sample
        curNeuronConcatMtx = curConcats(:, fannyShape).';
        curNeuronConcatMtx = reshape(curNeuronConcatMtx,[5,numTrials]);
        concatFannyShaped(fannyShape,:,:) = curNeuronConcatMtx; 
    end 

    concat_fanos = ifanofactor (concatFannyShaped);

    figure(201); 
    subplot(1,4,feedMeSeymour)
    boxplot([concat_fanos, concat_R_squared'])
    xticklabels({[curConcatName ' FFs'], [curConcatName ' r^2']})

    % correlation value
    concat_cor = corr(concat_fanos,concat_R_squared');
    
    % correlation figure 
    figure(202); 
    subplot(1,4,feedMeSeymour)
    title('Relationship between Fano Factor and R^2')
    
    scatter(concat_fanos,concat_R_squared)
    legend(['r = ' num2str(concat_cor)]);
    ylabel('R squared')
    xlabel('Fano Factor')
    title([curConcatName ' Neurons'])

    saveResults(feedMeSeymour).name = curConcatName;
    saveResults(feedMeSeymour).Rsquared = concat_R_squared;
    saveResults(feedMeSeymour).Predicted = Predicted; 
    saveResults(feedMeSeymour).fanoFactors = concat_fanos; 
    

end 

figure(200)
linkaxes
figure(201)
linkaxes
figure(202)
linkaxes



%% Multiple Regression Behavior from Neurons
selectedBehavior = trials.responseLatency; %response time per trial; use multiple neurons to predict 
RespSpikesPerTrial = squeeze(mean(selected_Responsive,2));

Predicted = iregressbehavior(spikes,selectedBehavior);

% correlation figure 
figure; 
title('Relationship between Fano Factor and R^2')

subplot(1,2,1)
scatter(resp_fanos,resp_R_squared)
legend(['r = ' num2str(resp_cor)]);
ylabel('R squared')
xlabel('Fano Factor')
subtitle('Reponsive Neurons')

%% Scaling up - apply to multiple brain areas 
region1 = 'VISp';
region2 = 'MOs';
region_index = find(strcmp(regions.name, region1));
region_neurons_1 = find(neurons.region == region_index);
region_index = find(strcmp(regions.name, region2));
region_neurons_2 = find(neurons.region == region_index);


[spike_counts_1, ~, allSpikesperTrial1] = get_spike_counts(trials, S, region_neurons_1); 
[spike_counts_2, ~, allSpikesperTrial2] = get_spike_counts(trials, S, region_neurons_2);

% pull out the r squared after regressing across two brain areas
[Predicted1, Predicted2,R_squared1,R_squared2, units_region1, units_region2] = ...
    TwoRegionRegress(spike_counts_1,spike_counts_2, num_sample);

%% Fano Factors across regions

%spike_counts_downsample = time x neurons 
% figures for responsive and unresponsive cell groups
fanos1 = ifanofactor (allSpikesperTrial1(units_region1,:,:));
fanos2 = ifanofactor (allSpikesperTrial2(units_region2,:,:));

% correlation value
fanos1 = fanos1(~isnan(fanos1));
R_squared1 = R_squared1(~isnan(R_squared1));
cor1 = corr(fanos1,R_squared1');
% remove nan spike 
fanos2 = fanos2(~isnan(fanos2));
R_squared2 = R_squared2(~isnan(R_squared2));
cor2 = corr(fanos2,R_squared2');

% correlation figure 
figure; 
title('Relationship between Fano Factor and R^2')

subplot(1,2,1)
scatter(fanos1,R_squared1)
legend(['r = ' num2str(cor1)]);
ylabel('R squared')
xlabel('Fano Factor')
subtitle('Region 1')

subplot(1,2,2)
scatter(fanos2,R_squared2)
legend(['r = ' num2str(cor2)]);
ylabel('R squared')
xlabel('Fano Factor')
subtitle('Region 2')

