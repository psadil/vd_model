function [selectivity, act_peak, act_out] = calc_selectivity(win_row, win_col, dist_mse,p)
% calc_selectivity -- calculates selectivity of grid to stimulus input.
% Currently, calculates both in terms of classic logistic but also now new
% guassian selectivity.

% called by: present_stimulus
% calls: NA

% input
%   win_row / win_col : indicies of winning row and winning col
%   dist_mat: distance of every node from stimulus input in mse space
%   p: experimental structure

% output:
%   selectivity: selectivity of this grid to the input
%   act_peak: activation of the nodes determined to be in the peak
%   act_total: total activation of entire grid

 

%%


% logistic based activation

% since dist_mat is in distance, it needs to be inverted to make smaller
% distances have greater activation. log is because it works (I think it
% works because it squashes the massive difference in distance from PRC and
% caudal layer)
act_out = log(ones(p.nRows,p.nRows) ./ dist_mse);     

% logistic function squashes activaiton to be no more than 1 (but also no
% less than .5)
act_out = 1./(1+exp(-p.k_expt*act_out));


%%

% find indicies of every potential winning node

% initialise array and record winner for all situations
winners = zeros(25,2);
winners(1,:) = [win_row win_col];
%city block distance 1 neighbours
winners(2,:) = [win_row win_col+1];
winners(3,:) = [win_row win_col-1];
winners(4,:) = [win_row+1 win_col];
winners(5,:) = [win_row-1 win_col];
%city block distance 2 neighbours
winners(6,:) = [win_row+1 win_col+1];
winners(7,:) = [win_row+1 win_col-1];
winners(8,:) = [win_row-1 win_col-1];
winners(9,:) = [win_row-1 win_col+1]; 
%city block distance 2-in-a-line neighbours
winners(10,:) = [win_row win_col+2]; 
winners(11,:) = [win_row win_col-2];
winners(12,:) = [win_row+2 win_col];
winners(13,:) = [win_row-2 win_col];
%city block distance 3 neighbours
winners(14,:) = [win_row+1 win_col+2]; 
winners(15,:) = [win_row+1 win_col-2];
winners(16,:) = [win_row-1 win_col+2];
winners(17,:) = [win_row-1 win_col-2];
winners(18,:) = [win_row+2 win_col+1]; 
winners(19,:) = [win_row+2 win_col-1];
winners(20,:) = [win_row-2 win_col+1];
winners(21,:) = [win_row-2 win_col-1];
%city block distance 3-in-a-line neighbours
winners(22,:) = [win_row win_col+3]; 
winners(23,:) = [win_row win_col-3];
winners(24,:) = [win_row+3 win_col];
winners(25,:) = [win_row-3 win_col];

for win_unit = 1:p.sizeOfPeak 
    if winners(win_unit,1) > p.nRows,
        winners(win_unit,1) = winners(win_unit,1) - p.nRows;
    elseif winners(win_unit,1) <= 0,
        winners(win_unit,1) = winners(win_unit,1) + p.nRows;
    end
    if winners(win_unit,2) > p.nRows,
        winners(win_unit,2) = winners(win_unit,2) - p.nRows;
    elseif winners(win_unit,2) <= 0,
        winners(win_unit,2) = winners(win_unit,2) + p.nRows;
    end
end


% calc activation of all nodes in actual peak
act_peak = 0;
for unit = 1:p.sizeOfPeak
    act_peak = act_peak + act_out(winners(unit,1), winners(unit,2));
end


selectivity = act_peak/(sum(act_out(:)));



end