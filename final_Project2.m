%% Directory path across computers

%% Load Helper Functions functions
code_dir = fullfile("C:", "Neuda2023", "Code", "w1d1", "npy-matlab-master");
addpath(genpath(code_dir))

%% Load data
data_dir = fullfile("C:", "Neuda2023", "Data", "Steinmetz", "allData");
cd(data_dir)
sesPath = "Moniz_2017-05-15";

% load all variables for session into a struct
S = loadSession(sesPath); % this calls a custom-coded read function for this dataset; it reads all .npy and .tsv files in the directory
clear ses sesName sesPath

%% list of known computers: point to git repo directory on your device
if exist('C:\Users\kedea\Documents\CSHL_GroupProject2','dir')
    home_dir = 'C:\Users\kedea\Documents\CSHL_GroupProject2';
elseif exist('C:\Users\course\Documents\CSHL_GroupProject2','dir')
    home_dir = 'C:\Users\course\Documents\CSHL_GroupProject2';
else
    error('Add your directory to this list (or rewrite this if there is a better way) - Kat')
% code_dir = fullfile("C:", "Neuda2023", "Code", "w1d1", "npy-matlab-master");
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
whichMouse = 'Hench_2017-06-18';

% genpath now sees all of the data inside allData
sesPath = [home_dir '\allData\' whichMouse];

% load all variables for session into a struct
S = loadSession(sesPath); % this calls a custom-coded read function for this dataset; it reads all .npy and .tsv files in the directory
clear ses sesName sesPath

%% Region Properties
% make a structure containing properties of the regions included in the session

regions = struct;
regions.name = unique(S.channels.brainLocation.allen_ontology,'rows'); % a character array 
regions.N = size(regions.name,1);
regions = orderfields(regions,{'N','name'});
regions.color = hsv(size(regions.name,1)-1); % unique rgb triplet for each region
regions.color(regions.N,:) = [ .5 .5 .5]; % grey for 'root' (i.e. not known)
% go further: add probe and depth fields to regions struct

%% Extract Neurons
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

[~,depthOrder ] = sort(neurons.depth);


%% Clean up data
%get rid of outliers



%% Select one area in animal



%% Multiple Regression Behavior from Neurons



%% Scaling up - apply to multiple brain areas 
