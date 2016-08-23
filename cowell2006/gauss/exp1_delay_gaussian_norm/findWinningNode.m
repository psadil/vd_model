function [ win_row, win_col, dist_mat ] = findWinningNode( weights, inp_mat )
%%% Find winning node, dot-product

% dist_mat = sqrt(sum((weights-inp_mat).^2,3));

% take dot product along weight dimension
dist_mat = dot(inp_mat,weights,3);

% Finds the row and column of max similarity
[win_row,win_col] = find(dist_mat==max(dist_mat(:)));


% If length(win_rows) is 2, chooses between first and second elements randomly
if length(win_row)>1 
    rand_idx = ceil(length(win_row)*rand);	
    win_row = win_row(rand_idx);
    win_col = win_col(rand_idx);
end


end

