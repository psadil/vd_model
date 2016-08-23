function [ win_row, win_col, dist_mat ] = findWinningNode( weights, inp_mat, nInpDims )
%%% Find winning node


%
% last set of working code below...

% OH, this is actually just the mean squared error, not euclidean distance at
% all!

dist_mat_slices_sq = (weights - inp_mat).^2; %% Use portion of stimulus corresponding to the current grid
% First slice of this contains (squared) distances of x-coords,
% second slice contains (squared) distances of y-coords
dist_mat_temp_sq = sum(dist_mat_slices_sq,3);

% Sum across x and y slices to get total distance
% dist_mat_slices = (weights - inp_mat);
% dist_mat_temp = sum(dist_mat_slices,3);

dist_mat = dist_mat_temp_sq/nInpDims;

% dist_mat_sqrt = sqrt(dist_mat_temp_sq);
% Normalise by dividing by num_input_dims
[win_row,win_col] = find(dist_mat==min(dist_mat(:)));
% Finds the row and column of minimal distance grid point(s)
if length(win_row)>1
    rand_idx = ceil(length(win_row)*rand);	% If length(win_rows) is 2, chooses between first and second elements randomly
    win_row = win_row(rand_idx);
    win_col = win_col(rand_idx);
end



end

