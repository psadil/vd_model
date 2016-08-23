function [p, stims] = createDelayStimuli(p)
% createDelayStimuli - creates stimuli for input to simulation of delay

% still want to make sure that no stimuli repeat

%%
% stims = zeros(max(p.nTrials), p.components,2);

% generate stims
stims = rand(max(p.nTrials), p.components,1);

% normalize
stims = stims ./ repmat(sqrt(sum(stims.^2,2)),[1,p.components]);

% make every second stim in the pair .5 away from the first
stims(:,:,2) = mod(stims(:,:,1)+p.shift,1);

% normalize
stims(:,:,2) = stims(:,:,2) ./ ...
    repmat(sqrt(sum(stims(:,:,2).^2,2)),[1,p.components]);


end