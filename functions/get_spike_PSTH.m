function spike_PSTH = get_spike_PSTH(trials, S, neurons_indices)

stimTimes = trials.visStimTime;
num_bin = length(-0.6:0.01:0.6);
num_trial = length(stimTimes);
num_unit = length(neurons_indices);


%% Count total spikes within regions at medium (100ms) time resolution - 2 sec
timeBins = nan(num_bin, length(stimTimes));

for i = 1 : length(stimTimes)
    timeBins(:, i) = double([(stimTimes(i) - 0.7) : 0.01 : (stimTimes(i) + 0.5)]);  
end

spike_PSTH = nan(num_unit, num_trial, num_bin-11);


% spike_counts  = zeros(length(timeEdges100ms)-1, length(neurons_indices));

for nNum = 1:length(neurons_indices) 
    
    nn = S.spikes.clusters == neurons_indices(nNum); % get indices of neurons in region
    this_neuron = S.spikes.times(nn);
    
    for ntrial = 1 : num_trial
        aaa = conv(histcounts(this_neuron, timeBins(:, ntrial)), ones(1,10))';
        spike_PSTH(nNum, ntrial, :) = aaa(11:end-9);

    end
    
end


end