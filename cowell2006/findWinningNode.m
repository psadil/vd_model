function [ win_row, win_col, dist_mse ] = findWinningNode( weights, inp )
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
dist_mat_sse = sum((weights - inp).^2,3);

% mean of summed squared distance
dist_mse = dist_mat_sse/size(weights,3);

% BUT, it's euclidean distance that's used to find winning node.
dist_euclid = sqrt(dist_mat_sse);

% find that winning node
[win_row,win_col] = find(dist_euclid==min(dist_euclid(:)));

% incase two winning nodes were found, pick only one of them. Rarely, if
% ever, happens
if length(win_row) > 1
    rand_idx = ceil(length(win_row)*rand);
    win_row = win_row(rand_idx);
    win_col = win_col(rand_idx);
end

end

