function [ win_row, win_col, dist_mat ] = findWinningNode( weights, inp_mat, nInpDims )
%%% Find winning node
% this function calculates two things. first is the winning node, with is
% the node that is the closest (in a euclidean sense) to the input.

% ALSO, the function calculates dist_mat, which is the distance of the rest
% of the grid's nodes distance from the winning node (distance in the
% city-block sense).


dist_mat_slices_sq = (weights - inp_mat).^2; %% Use portion of stimulus corresponding to the current grid

% First slice of this contains (squared) distances of x-coords,
% second slice contains (squared) distances of y-coords
dist_mat_temp_sq = sum(dist_mat_slices_sq,3);

% Sum across x and y slices to get total distance
% dist_mat_slices = (weights - inp_mat);
% dist_mat_temp = sum(dist_mat_slices,3);
 
dist_mat = dist_mat_temp_sq/nInpDims;

dist_mat_sqrt = sqrt(dist_mat_temp_sq);
% Normalise by dividing by num_input_dims
[win_rows,win_cols] = find(dist_mat_sqrt==min(min(dist_mat_sqrt)));
% Finds the row and column of minimal distance grid point(s)
rand_idx = ceil(length(win_rows)*rand);	% If length(win_rows) is 2, chooses between first and second elements randomly
win_row = win_rows(rand_idx);
win_col = win_cols(rand_idx);

end

