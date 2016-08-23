function [ predictions ] = getDistance( parms )
%getDistance calculates average first and second half distances
%   I think this function is for a simplier version of the model, one that
%   doesn't apply to the simplexGO folder

global consts;


dist0 = abs(consts.input - consts.w0);

% estimate halflife parameter or decay constant

v_r = exp(-((consts.r ./ parms(2)).^2)); % should just be 1
f_r = parms(1) * v_r;

dist0_full = mean(dist0.^2,2);

act0 = log(1./dist0_full);
act0_squashed = 1./(1+exp(-1*act0));
act0_squashed_from1 = 1-act0_squashed;

weights = zeros(consts.nTrials,25,3);
weights(1,:,:) = consts.w0;

diffFromInp = zeros(consts.nTrials,25,3);
diffFromInp(1,:,:) = dist0;



trial = 1;
ratio = 1;
while mean(ratio) > .5
    trial = trial + 1;
    
    diffFromInp(trial,:,:) = consts.input - squeeze(weights(trial-1,:,:));
    weights(trial,:,:) = squeeze(weights(trial-1,:,:)) + repmat(f_r',[1,3]).*(squeeze(diffFromInp(trial,:,:)));
    
    dist_new = mean(squeeze(diffFromInp(trial,:,:).^2),2);
    act_new = log(1 ./ dist_new);
    
    act_squashed = 1./(1+exp(-1*act_new));
    act_squashed_from1 = 1-act_squashed;
    ratio = mean(act_squashed_from1) ./ mean(act0_squashed_from1);
    
end

tau_trial = trial-2;

% half_d0 = d0/2;
%
% [~, tau_trial] = min(abs(half_d0 - diffFromInp(2:end)));
% % tau_trial = tau_trial;
%
% tau_trial = 1;
%% distances
trials = 1:consts.nTrials;

% need to pull out activations now that dealing with vectors!!!

dist_LA = act0_squashed * (1/2).^((consts.trialLA * trials)/tau_trial);
dist_HA = act0_squashed * (1/2).^((consts.trialHA * trials)/tau_trial);


% distance of averages, LA
dist_LA_first = mean(dist_LA(1:36));
dist_LA_second = mean(dist_LA(37:72));

diff_LA = dist_LA_first-dist_LA_second;

% distance of averages, HA
dist_HA_first = mean(dist_HA(1:36));
dist_HA_second = mean(dist_HA(37:72));

diff_HA = dist_HA_first-dist_HA_second;

predictions = [diff_LA, diff_HA];

end

