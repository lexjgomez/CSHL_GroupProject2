%% pulll out responsive cells

function [responsiveUnits, nonresponsiveUnits, Indices] = classify_by_response(allSpikesperTrial,num_unit_all)

preStim = allSpikesperTrial(:,1:5,:);
postStim = allSpikesperTrial(:,6:10,:);

%do t-test
for th = 1:num_unit_all
    [h,p] = ttest2(mean(preStim(th,:,:),3),mean(postStim(th,1:5,:),3));
    storeTests(th) = p;
    if p < 0.05
       checkTail = mean(mean(preStim(th,:,:),3))<mean(mean(postStim(th,:,:),3)); %checks to see if the spikec count during stim time was larger than the baseline
       if checkTail == 1
        check_all_Tails(th) = 1;
       else
        check_all_Tails(th) = 0;
       end 
    else
        check_all_Tails(th) = 0;
    end 
end 

Indices.Trial_IndiciesResponsive = find(check_all_Tails);
    check_all_Tails_inv = abs(check_all_Tails-1);
Indices.Trial_IndiciesNonResponsive = find(check_all_Tails_inv);



%pull out cells whose fr significantly increased after stimulus
responsiveUnits = allSpikesperTrial((storeTests<0.05 & check_all_Tails>0),:,:);
negResponse = allSpikesperTrial((storeTests<0.05 & check_all_Tails<0),:,:);
%basically everything else
nonresponsiveUnits = allSpikesperTrial((storeTests>0.05),:,:);
