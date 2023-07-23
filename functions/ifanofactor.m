function fanos = ifanofactor(allSpikes_downsample)
%fano factors that are above 1 means that you are getting more variability,
%and getting different number of spikes per trial

%low fano factor means less variability so you should get the same ish
%number of spikes per trial, and if this occurs during the stimulus, that
%makes sense because the neuron is responding to the stimulus #synchronized
% :)

for ii = 1:size(allSpikes_downsample,1) %Going through all relevant neurons
   currNeuronSpikes = allSpikes_downsample(ii,:);
   neuronMeanVar(ii,1) = mean(currNeuronSpikes); %Trial mean
   neuronMeanVar(ii,2) = var(currNeuronSpikes); %Trial by trial variability
   fanos(ii,1) = neuronMeanVar(ii,2)/neuronMeanVar(ii,1); %Fano factor
end

%Plotting it
figure
begRange = 1e-3;
endRange = 1e3;
h = plot(neuronMeanVar(:,1),neuronMeanVar(:,2),'.','MarkerSize',20,'color','k');
line([begRange endRange],[begRange endRange],'color','b'); %Poisson expectation
xlim([begRange endRange]); ylim([begRange endRange]); 
xlabel('Spikecount mean')
ylabel('Spikecount variance')
movshonize(26,1)
set(gca,'XScale','log')
set(gca,'YScale','log')

%Now make a histogram of the fano factors in this region
%Not really that meaningful, as what matters is the proportion above and
%below 1
figure
nBins = 31; 
h = histogram(fanos,nBins);
h.FaceColor = 'k';
h.FaceAlpha = 1;
h.EdgeColor = 'w';
ylim([0 max(h.Values)+1])
line([1 1], [0 max(h.Values)+1],'color', 'r', 'linewidth', 3 )
xlabel('Fano factor')
ylabel('Number of neurons per bin')
movshonize(26,1); makeWhite
