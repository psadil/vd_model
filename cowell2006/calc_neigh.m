function f_out = calc_neigh(win_row, win_col, layer, p, delay)
% calc_neigh --  determinings the neighborhood.
% This defines the spire around which grids will be updated.

% called by: pretrain, present_stimulus
% calls: NA

% input:
%   win_row / win_col -- indices of winning node
%   dist_mat: distance of every node away from sample stimulus in mse
%      space.
%   layer: layer in which the current grid comes from
%   p: experimental structure

% output
%   f_out: neighborhood function. defines by how much nodes will be updated
%      to more closely match stimulus input. numRows x numRows x
%      numInputDims(layer)
%   act_out: not currently utilized
% Calculate city-block distance from winner in grid, with wraparound

%%
% find distance of each unit from winner (using grid_matrix, which stores the position of each unit in the grid)
% create a matrix with a slice for each of the two potential minimum distances for rows and cols

%% grab city-block distance of every node away from winning node

% first, look at row and column distance. Because the grid 'wraps around'
% into a toroid, there could be two possible lowest distances in each dim
row_dist_mat(:,:,1) = abs(p.gridMat(:,:,1) - win_row);
row_dist_mat(:,:,2) = p.nRows - abs(p.gridMat(:,:,1) - win_row);
col_dist_mat(:,:,1) = abs(p.gridMat(:,:,2) - win_col);
col_dist_mat(:,:,2) = p.nRows - abs(p.gridMat(:,:,2) - win_col);

% find the minimum of the two possible distances for each row and col
min_row_dist_mat = min(row_dist_mat,[],3);
min_col_dist_mat = min(col_dist_mat,[],3);

% Sum the two minimum distances to get overall city_block distance
grid_dist = min_row_dist_mat + min_col_dist_mat;


%% Calculate Gaussian movement-strength function for each node

if delay
    eta = p.eta_int;
else
    eta = p.eta;
end


f_1dim = eta*exp(-(grid_dist/p.G).^2);

% need a gaussian for each dimension in attended to by the grid.
f_out = repmat(f_1dim, [1,1,p.nInputDims(layer)]);


end
