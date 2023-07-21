function spike_counts = get_spike_counts(trials, S, neurons_indices)
stimTimes = trials.visStimTime; 


%% Count total spikes within regions at medium (100ms) time resolution - 2 sec
timeEdges100ms = 0:.1:(stimTimes(end)+3) ; % end points of 100 ms bins up to end of trials period
spike_counts  = zeros(length(timeEdges100ms)-1, length(neurons_indices));
for nNum = 1:length(neurons_indices) % last region is usually 'root' (i.e. unassigned)
    nn = S.spikes.clusters == neurons_indices(nNum); % get indices of neurons in region 
    spike_counts(:,nNum) = histcounts(S.spikes.times(nn),timeEdges100ms)'; % complete set of bins
end
