function [Predicted1, Predicted2,R_squared1,R_squared2, units_region1,units_region2] = ...
    TwoRegionRegress(spike_counts_1,spike_counts_2, num_sample)


sampleNum_1 = size(spike_counts_1,2);
units_region1 = randsample(1:sampleNum_1, num_sample);
spike_counts_1 = spike_counts_1(:, units_region1);
Loss_1 = nan(1, num_sample);
Predicted1 = nan(size(spike_counts_1));


sampleNum_2 = size(spike_counts_2,2);
units_region2 = randsample(1:sampleNum_2, num_sample);
spike_counts_2 = spike_counts_2(:, units_region2);
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


RSS = (spike_counts_1 - Predicted1) .^ 2;
RSS = mean(RSS, 1);
TSS = (spike_counts_1 - mean(spike_counts_1, 1)) .^ 2;
TSS = mean(TSS, 1);

R_squared1 = 1 - (RSS./TSS);

figure;
subplot(1,2,1);
histogram(R_squared1);
xlabel('R squared');
ylabel('Number of Neurons');
title('Brain Region 1 (Predicted by region 2)');
set(gca, 'box', 'off');
set(gca, 'tickdir', 'out');



RSS = (spike_counts_2 - Predicted2) .^ 2;
RSS = mean(RSS, 1);
TSS = (spike_counts_2 - mean(spike_counts_2, 1)) .^ 2;
TSS = mean(TSS, 1);

R_squared2 = 1 - (RSS./TSS);

subplot(1,2,2);
histogram(R_squared2);
xlabel('R squared');
ylabel('Number of Neurons');
title('Brain Region 2 (Predicted by region 1)');
set(gca, 'box', 'off');
set(gca, 'tickdir', 'out');


end
