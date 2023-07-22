function Predicted = imultipleregress(spike_counts_downsample,num_sample)


Loss = nan(1, num_sample);
Predicted = nan(size(spike_counts_downsample));
for i = 1 : num_sample

    Y = spike_counts_downsample(:, i);
    X = spike_counts_downsample;
    X(:, i) = [];

    % cross validation with lasso regularization
    Mdl = fitrlinear(X,Y, 'Learner','leastsquares','CrossVal','on','Regularization','lasso');
    Loss(i) = kfoldLoss(Mdl);
    Predicted(:, i) = kfoldPredict(Mdl);
end

