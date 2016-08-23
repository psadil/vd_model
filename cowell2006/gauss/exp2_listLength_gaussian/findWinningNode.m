function [ win_row, win_col, dist_mat ] = findWinningNode( weights, inp_mat )
% FindWinningNode -- finds node in available grid that best matches current
% input.
% NOTE: this would be a place to adjust if continuing to use SOMs. Current
%    best node is determined by minimum euclidean distance. Alternatives
%    include highest dot product.


% called by: pretrain, present_simulus
% calls: NA

% input:
%    weights: weights of current grid to which the input will be compared
%       against.
%   inp_mat: stimulus to which the weights will be compared against.
% NOTE: inputs should be of same size

% output:
%   win_row / win_col: indices of winning node in grid
%   dist_mat: distance of every node in grid away from stimulus input.
%      Measured in mean-squared-error, NOT euclidean distance. This will be
%      used in calculations of activation and selectivity in
%      calc_selectivity


%%

% squared distance of inputs from weights
dist_mat_slices_sq = (weights - inp_mat).^2;

% sum of squared distance
dist_mat_temp_sq = sum(dist_mat_slices_sq,3);

% mean of summed squared distance
dist_mat = dist_mat_temp_sq/size(weights,3);

% BUT, it's euclidean distance that's used to find winning node.
dist_mat_sqrt = sqrt(dist_mat_temp_sq);

% find that winning node
[win_row,win_col] = find(dist_mat_sqrt==min(min(dist_mat_sqrt)));

% incase two winning nodes were found, pick only one of them. Rarely, if
% ever, happens
if length(win_row) > 1
    rand_idx = ceil(length(win_row)*rand);
    win_row = win_row(rand_idx);
    win_col = win_cols(rand_idx);
end

end

