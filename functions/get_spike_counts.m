function [spike_counts, allSpikes, allSpikesperTrial] = get_spike_counts(trials, S, neurons_indices)
stimTimes = trials.visStimTime; 


%% Count total spikes within regions at medium (100ms) time resolution - 2 sec
timeEdges100ms = 0:.1:(stimTimes(end)+3) ; % end points of 100 ms bins up to end of trials period
spike_counts  = zeros(length(timeEdges100ms)-1, length(neurons_indices));
win_before = 5;  %0.5seconds 
win_after = 5; %0.5 seconds
allSpikes = []; %array for all spikes per trial per neuron
allSpikesperTrial = [];

for nNum = 1:length(neurons_indices) % last region is usually 'root' (i.e. unassigned)
    %get spikes for every 0.1ms timebin
    nn = S.spikes.clusters == neurons_indices(nNum); % get indices of neurons in region 
    spike_counts(:,nNum) = histcounts(S.spikes.times(nn),timeEdges100ms)'; % complete set of bins

    temp = {};
    % get the total spikes per trial 
    for i = 1:size(stimTimes,1)
        currStim = stimTimes(i);
        [~, idx] = min(abs(timeEdges100ms-currStim));
        idxBefore = idx - win_before; idxAfter = idx + win_after;
        currSpikesNeuron = spike_counts(:,nNum);

        %neuron x time x trial
        allSpikesperTrial(nNum, :,i) = currSpikesNeuron(idxBefore:idxAfter)';

        %neuron spikes x trial
        allSpikes(nNum,i) = sum(currSpikesNeuron(idxBefore:idxAfter));
    end 

end


