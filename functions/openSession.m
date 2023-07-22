function [S, neurons, trials, regions] = openSession(sesPath)
%% Specify the path to the folder containing the datasets
if ~exist('sesPath'), sesPath = 'Moniz-05-18'; end % sample with both motor and sensory areas

% specify path to data

% load all variables for session into a struct
S = loadSession(sesPath); % this calls a custom-coded read function for this dataset; it reads all .npy and .tsv files in the directory
% LFP .bin files and .csv files from DLC must be read in separately
% NB sometimes directory is missing a few .npy files; also a few are mixed up (e.g. Moniz-5-18 spike pattern length)
clear ses sesName sesPath

%% Could load several sessions
% mainPath = './'; % CHANGE 
% % generate list of available sessions
% sessionList = dir(mainPath);
% sessionList = sessionList(4:end);  % For Mac: first 3 entries are directories or hidden file

%% Region Properties
% make a structure containing properties of the regions included in the session

regions = struct;
regions.name = unique(S.channels.brainLocation.allen_ontology,'rows'); % a character array 
regions.N = size(regions.name,1);
regions = orderfields(regions,{'N','name'});
regions.color = hsv(size(regions.name,1)-1); % unique rgb triplet for each region
regions.color(regions.N,:) = [ .5 .5 .5]; % grey for 'root' (i.e. not known)
% go further: add probe and depth fields to regions struct

%% Neuron Properties 
% neuron ('cluster') properties are distributed in a few different places,
% we bring them together into a new struct called 'neurons'

neurons = struct;
% extract the neuron properties
neurons.id = unique(S.spikes.clusters);
% note that id's in S.spikes.clusters run from 0 to (N-1), and do not match row index of S.clusters 
neurons.N = length(neurons.id); % number of neurons
neurons = orderfields(neurons,{'N','id'});

% identify region by row in region struct
[~,Loc] = ismember(S.channels.brainLocation.allen_ontology(S.clusters.peakChannel,:),regions.name,'rows');
neurons.region = Loc; % numeric code for region
regions.name = strtrim(string(regions.name)); % after match, make names into a more manageable string variable
clear Loc
neurons.depth = S.clusters.depths;
neurons.probe = S.clusters.probes;


%% Create a structure containing information about the trial
trials = struct;
trials.N = size(S.trials.intervals,1);
trials.isStimulus = S.trials.visualStim_contrastLeft > 0 | S.trials.visualStim_contrastRight > 0; % did a brighter image appear on one side?
trials.isMovement = S.trials.response_choice ~= 0; % did the mouse move the wheel in response?
trials.visStimTime = S.trials.visualStim_times; % time of stimulus
trials.responseTime = S.trials.response_times; % time response recorded
trials.turn = S.trials.response_choice; % 
trials.contrast =  S.trials.visualStim_contrastLeft - S.trials.visualStim_contrastRight ; % contr
LvsR=sign( trials.contrast ); % indicates contrast difference relevant for choice
xtab1=crosstab( LvsR ,  S.trials.response_choice);
% stimuli  L<R, L=R, L>R in rows; choices -1,0,1 in column
% Make a contingency table of stimuli x choices and display it on console
table(xtab1(:,1),xtab1(:,2),xtab1(:,3),'VariableNames',{ 'turn R', 'No turn', 'turn L'},'RowNames',{'L', 'equal', 'R'});
trials.Correct = ( LvsR == 0 | S.trials.response_choice == LvsR) & trials.isMovement ; % flag correct choices
% omitting passive if both 0
trials.responseLatency = trials.responseTime - S.trials.goCue_times; 
trials.timeOut = trials.responseLatency > 1.49; % time outs at 1.5 sec
