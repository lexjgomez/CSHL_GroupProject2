function iregressbehavior(allSpikes_downsample,selectedBehavior)

predictors = allSpikes_downsample;
outcome = selectedBehavior;
numTrials = size(allSpikes_downsample,2);
numNeurons = size(allSpikes_downsample,1);

Pred = [ones(numTrials,1) predictors']; % the predictors ; add ones for intercept

%visualize the R2 value first
%lasso(outcome, Pred)
% [b, bint, r, rint,stats] = regress(outcome, Pred);
lasso()
figure()
bar(b(2:length(b)))
xlabel('Neurons'); ylabel('b');
title('Multiple Regression: Neuronal Activity vs Response Latency')
subtitle(strcat("R2 = ", num2str(stats(1))))

%create a linear regression model
Mdl = fitrlinear(Pred, outcome, 'Learner','leastsquares','CrossVal','on','Regularization','lasso');
Loss = kfoldLoss(Mdl); %because this is a fitrlinear model, 1-R2, with adjustment term.  %loss should be 0 - 1, 1 = predict nothing
Predicted = kfoldPredict(Mdl); %predict the results, predict the training data, 

RSS = mean(outcome - Predicted .^ 2,1);
TSS = mean((outcome - mean(outcome, 1)) .^ 2,1);
R_squared = 1 - (RSS./TSS);

% % % doing regression on random values
% % Pred2 =  [ ones(numTrials,1)  randn(numTrials,numNeurons)];
% % [b, bint, r, rint,stats] = regress(outcome, Pred2);
% % stats(1) % R^2 
% % figure()
% % bar(b(2:length(b)))
% % xlabel('Neurons')
% % ylabel('b');
% % title('Multiple Regression: Random Neurons vs Response Latency')




%% Ways to reduce over-fitting
% default LASSO to prune down likely list of neurons
% Pred = [ones(numTrials,1) predictors']; % the predictors ; add ones for intercept
% B = lasso( Pred, outcome' );
% % B contains coefficients that are estimated with various penalties (default 100)
% figure
% numNonZeros = sum(B~=0); % adds up all zero coefficients
% plot(numNonZeros);  % where do you get ~50 non-zero coefficients 
% yline(50)
% figure; bar(B(:,90))
% 
% % Let's look inside
% [B, fitInfo] = lasso( Pred, data.Eye_full(1,:)' );
% % What values of lambda the reg parameter did MATLAB use?
% plot(fitInfo.Lambda) % looks non-linear
% plot(log(fitInfo.Lambda)) % MATLAB picks parameters evenly by logarithm
% plot(fitInfo.MSE) % goes from


% where do you draw the line?
% What's a reasonable error?
