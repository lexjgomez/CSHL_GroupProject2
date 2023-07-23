function [Predicted1, Predicted2] = TwoRegionRegress(spike_counts_1,spike_counts_2, num_sample)


sampleNum_1 = size(spike_counts_1,2);
spike_counts_1 = spike_counts_1(:, randsample(1:sampleNum_1, num_sample));
Loss_1 = nan(1, num_sample);
Predicted1 = nan(size(spike_counts_1));


sampleNum_2 = size(spike_counts_2,2);
spike_counts_2 = spike_counts_2(:, randsample(1:sampleNum_2, num_sample));
Loss_2 = nan(1, num_sample);
Predicted2 = nan(size(spike_counts_2));




for i = 1 : num_sample

    Y = spike_counts_1(:, i);
    X = spike_counts_2;

    % cross validation with lasso regularization
    Mdl_1 = fitrlinear(X,Y, 'Learner','leastsquares','CrossVal','on','Regularization','lasso');
    Loss_1(i) = kfoldLoss(Mdl_1);
    Predicted1(:, i) = kfoldPredict(Mdl_1);
    
    
    Y = spike_counts_2(:, i);
    X = spike_counts_1;
    
    Mdl_2 = fitrlinear(X,Y, 'Learner','leastsquares','CrossVal','on','Regularization','lasso');
    Loss_2(i) = kfoldLoss(Mdl_2);
    Predicted2(:, i) = kfoldPredict(Mdl_2);
end






end
