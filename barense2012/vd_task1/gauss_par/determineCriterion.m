function [ p ] = determineCriterion( p, trial, threshUpdater )
%update criterion for use on next trial

% when two layers are available, each are allowed seperate decision
% thresholds

p.famil_diff_thresh = [threshUpdater , p.famil_diff_thresh(1:p.layer,1:end-1)];

end

