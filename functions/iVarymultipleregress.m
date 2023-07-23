function R_squared = iVarymultipleregress(spike_counts_downsample)

num_sample = size(spike_counts_downsample,2);
R_squared = nan(num_sample,num_sample-1);

for n= 1:num_sample %loop through all the neurons
        disp(num2str(n))
        specialNeuron = n;
        otherNeurons = setdiff([1:num_sample], specialNeuron); 
        parfor i = 1 : num_sample-1 %loop for each neuron, vary the number of predictors
            useNeurons = otherNeurons(1:i); %get other neurons in random sizes 
            outcome = spike_counts_downsample(:, specialNeuron);
            predictor = spike_counts_downsample(:,useNeurons);
            % cross validation with lasso regularization
            [Mdl] = fitrlinear(predictor, outcome, 'Learner','leastsquares','CrossVal','on','Regularization','lasso');
            Predicted = kfoldPredict(Mdl);
            
            RSS = mean((outcome - Predicted) .^ 2);
            TSS = mean((outcome - mean(outcome, 1)) .^ 2);
            R_squared(n,i) = 1 - (RSS./TSS);
        end
end 




