function [f_out, act_out] = calc_act_fast(win_row, win_col, dist_mat, layer, p,interefere)

% This function is called from delay_present_stimulus.m. Calculates the activation on the grid and the f matrix determining
% amount of learning that will accrue to each unit.
% Calculate city-block distance from winner in grid, with wraparound
if interefere
   alpha = p.eta_int; 
else
    alpha = p.eta;
end

%% find distance of each unit from winner (using grid_matrix, which stores the position of each unit in the grid) 
% create a matrix with a slice for each of the two potential minimum distances for rows and cols
row_dist_mat(:,:,1) = abs(p.gridMat(:,:,1) - win_row);
row_dist_mat(:,:,2) = p.numRows - abs(p.gridMat(:,:,1) - win_row);
col_dist_mat(:,:,1) = abs(p.gridMat(:,:,2) - win_col);
col_dist_mat(:,:,2) = p.numRows - abs(p.gridMat(:,:,2) - win_col);
%%find the minimum of the two values for row and for col
min_row_dist_mat = min(row_dist_mat,[],3);
min_col_dist_mat = min(col_dist_mat,[],3);
% Sum the two minimum distances to get the city_block distance
grid_dist = min_row_dist_mat + min_col_dist_mat;

% Calculate Gaussian movement-strength function for each node
f_1dim = alpha*exp(-(grid_dist/p.G).^2);
% f_1dim = f_1dim .* (grid_dist < p.filtPeak);

nInpDims=p.numInputDims(layer);
f_out=zeros(p.numRows,p.numRows,nInpDims);
for dim=1:nInpDims
    f_out(:,:,dim) = f_1dim;
end


act_out = log(ones(p.numRows,p.numRows)./dist_mat);
act_out = 1./(1+exp(-p.k_expt*act_out)); %squashing function

% act_out(act_out>9.21)=9.21;
 

end
     