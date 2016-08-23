function [ win_row, win_col, dist_mat ] = findWinningNode( weights, inp_mat, nInpDims )




dist_mat_temp_sq = sum((weights - inp_mat).^2,3); %% Use portion of stimulus corresponding to the current grid
% First slice of this contains (squared) distances of x-coords,
% second slice contains (squared) distances of y-coords
 
dist_mat = dist_mat_temp_sq/nInpDims;

dist_mat_sqrt = sqrt(dist_mat_temp_sq);
% Normalise by dividing by num_input_dims
[win_rows,win_cols] = find(dist_mat_sqrt==min(dist_mat_sqrt(:)));
% Finds the row and column of minimal distance grid point(s)
rand_idx = ceil(length(win_rows)*rand);	% If length(win_rows) is 2, chooses between first and second elements randomly
win_row = win_rows(rand_idx);
win_col = win_cols(rand_idx);


end