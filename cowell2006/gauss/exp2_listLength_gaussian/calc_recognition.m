function [ p ] = calc_recognition( p, selec_forComp, trial,gauss,stimSet)
% calc_recognition: calcs recognition socres
% NOTE: recognition is defined by equations 14-16 of cowell et al. 2006

% called by: model
% calls: NA

% input:
%   p: experimental structure
%   selec_forComp: selectivity from which to calculate recognnition score
%   trial: current trial (choice phase)
%   gauss: flag to indicate whether selec_forComp indicates logistic or
%      gaussian (when == 1) selectivity
%   stimSet: which stimulus set are we working on?

% output
%   p: now contains recognition score for this trial in this stimulus set.


%%

% if handling the control network, s_samp and s_nov are calculated as
% average across all available layers.
if p.layer == 2
    
    if gauss
        tmp_prc = squeeze(selec_forComp(2,1,:));
        tmp_caudal = squeeze(mean(selec_forComp(1,:,:),2));
        p.recogByLayer_gauss(trial,2,stimSet) = (tmp_prc(1) - tmp_prc(2)) / (tmp_prc(1) + tmp_prc(2));
        selec = mean([tmp_prc,tmp_caudal],2);
    else
        tmp_prc = squeeze(selec_forComp(2,1,:));
        tmp_caudal = squeeze(mean(selec_forComp(1,:,:),2));
        p.recogByLayer(trial,2,stimSet) = (tmp_prc(1) - tmp_prc(2)) / (tmp_prc(1) + tmp_prc(2));
        selec = mean([tmp_prc,tmp_caudal],2);
    end
else
    selec = squeeze(mean(selec_forComp(1,:,:),2));
end

if gauss
    tmp_caudal = squeeze(mean(selec_forComp(1,:,:),2));
    p.recogByLayer_gauss(trial,1,stimSet) = (tmp_caudal(1) - tmp_caudal(2)) / (tmp_caudal(1) + tmp_caudal(2));
    p.recognition_gauss(trial,stimSet) = (selec(1) - selec(2)) / (selec(1) + selec(2));
else
    tmp_caudal = squeeze(mean(selec_forComp(1,:,:),2));
    p.recogByLayer(trial,1,stimSet) = (tmp_caudal(1) - tmp_caudal(2)) / (tmp_caudal(1) + tmp_caudal(2));
    p.recognition(trial,stimSet) = (selec(1) - selec(2)) / (selec(1) + selec(2));
end

end