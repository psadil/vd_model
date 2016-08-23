function [ p ] = determineCriterion( p, trial )
%update criterion for use on next trial

% when two layers are available, each are allowed seperate decision
% thresholds

% p.mean_famil_diff(trial) = mean(familDiffs);
if p.numThresh == 2
    if all(p.usePRC(:,trial))
        famil_diff_thresh_new = p.famil_diff_thresh(2,1:end-1);
        p.famil_diff_thresh = [p.famil_diff_thresh(1,:); p.famil_difference(trial), famil_diff_thresh_new];
    else
        famil_diff_thresh_new = p.famil_diff_thresh(1,1:end-1);
        p.famil_diff_thresh = [p.famil_difference(trial), famil_diff_thresh_new; p.famil_diff_thresh(2,:)];
        % p.famil_diff_thresh = [p.famil_diff_thresh, p.famil_difference(trial)];
        %     p.threshForPlotting(trial) = mean(p.famil_diff_thresh);
        % p.familDiff_threshTracking(trial) = famil_diff_thresh_new;
    end
else
    famil_diff_thresh_new = p.famil_diff_thresh(1,1:end-1);
    p.famil_diff_thresh = [p.famil_difference(trial), famil_diff_thresh_new; p.famil_diff_thresh(2,:)];
    
end

end

