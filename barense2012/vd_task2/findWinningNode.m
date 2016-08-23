function [ win_row, win_col, dist_mat ] = findWinningNode( weights, inp_mat, nInpDims )
%%% Find winning node


% dist_mat_slices = (weights - inp_mat).^2;  % should not even be squared, according to Rosie's paper
% % % First slice of this contains (squared) distances of x-coords,
% % % second slice contains (squared) distances of y-coords
% dist_mat_temp = sum(dist_mat_slices,3);
% % Sum across the 2 matrix slices (or slices of all dimensions) to get total distance
%
% 
%
%
% dist_mat = dist_mat_temp./nInpDims;
% % Normalise by dividing by num_input_dims
%
% % dist_mat = sqrt(dist_mat_temp);
%
% % NEW CODE: finish distance calc with sqrt
% % dist_mat = sqrt(dist_mat);
%
% % [win_rows,win_cols] = find(dist_mat==min(min(dist_mat)));
% % [~, index] = min(sqrt(dist_mat_temp(:)));
% % [win_row, win_col] = ind2sub(size(dist_mat_temp),index);
% dist_mat_sqrt = sqrt(dist_mat_temp);
%
% [~,win_row] = min(min(dist_mat_sqrt,[],1));
% [~,win_col] = min(min(dist_mat_sqrt,[],2));
%
% if length(win_row)>1
%     % Finds the row and column of minimal distance grid point(s)
%     rand_idx = ceil(length(win_row)*rand);	% If length(win_rows) is 2, chooses
%     % between first and second elements randomly
%     win_row = win_row(rand_idx);
%     win_col = win_col(rand_idx);
% end
% %   if CYCLE == NUM_TRAINING_CYCLES,
% %      fprintf('\nVD_pretrain.m');
% %     fprintf('\ndist_mat');
% %    max(max(dist_mat))
% %   min(min(dist_mat))
% % end

%%
% last set of working code below...

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


%%
% % current set of working code below...
% 
% dist_mat_slices = (weights - inp_mat).^2; %% Use portion of stimulus corresponding to the current grid
% % First slice of this contains (squared) distances of x-coords,
% % second slice contains (squared) distances of y-coords
% dist_mat_temp = sum(dist_mat_slices,3);
% 
% % Sum across x and y slices to get total distance
% 
% 
% dist_mat_sqrt = sqrt(dist_mat_temp);
% 
% dist_mat = dist_mat_sqrt;
% 
% % Normalise by dividing by num_input_dims
% [win_rows,win_cols] = find(dist_mat_sqrt==min(min(dist_mat_sqrt)));
% % Finds the row and column of minimal distance grid point(s)
% rand_idx = ceil(length(win_rows)*rand);	% If length(win_rows) is 2, chooses between first and second elements randomly
% win_row = win_rows(rand_idx);
% win_col = win_cols(rand_idx);


end

