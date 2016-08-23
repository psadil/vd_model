function [ p ] = calc_corr( p, acts_forComp, trial, stimSet)
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
    
    tmp_prc = squeeze(acts_forComp(2,:,:,:,:));
    first_prc = tmp_prc(:,:,:,1);
    second_prc = tmp_prc(:,:,:,2);
    
    % calculate correlation
    rho_prc = corr(first_prc(:),second_prc(:), 'type','Pearson');
    
    % fisher transform
    p.corrByLayer(trial,2,stimSet) = .5*log((1+rho_prc)/(1-rho_prc));
    
    
end


% always calculate correlation for caudal layer
tmp_caudal = squeeze(acts_forComp(1,:,:,:,:));
first_caudal = tmp_caudal(:,:,:,1);
second_caudal = tmp_caudal(:,:,:,2);

rho_caudal = corr(first_caudal(:),second_caudal(:), 'type','Pearson');
p.corrByLayer(trial,1,stimSet) = .5*log((1+rho_caudal)/(1-rho_caudal));


% calculate correlation across all available layers
first = acts_forComp(1:p.layer,:,:,:,1);
second = acts_forComp(1:p.layer,:,:,:,2);

rho = corr(first(:),second(:), 'type','Pearson');
p.corr(trial,stimSet) = .5*log((1+rho)/(1-rho));


end