function [W, initial_selec, p] = present_stimulus(stim, W, p, trial)

inp_mat = repmat(reshape(stim,[1 1 length(stim)]), [p.nRows p.nRows 1]);

initial_selec = zeros(p.numLayers,max(p.numGrids));



%% Expose network to stimuli and update weights

for layer=1:p.layer
    
    firstFeatureToCheck=(1:p.numInputDims(layer):p.nGrids(layer)*p.numInputDims(layer));
    lastFeatureToCheck = firstFeatureToCheck+p.numInputDims(layer)-1;
    
    for grid = 1:p.nGrids(layer),
        

        input_mat=inp_mat(:,:,(firstFeatureToCheck(grid):lastFeatureToCheck(grid)));
        
        % put the variable 'weights' into the format previously accepted by the model.
        weights = squeeze(W(layer,:,:,1:p.numInputDims(layer),grid));
        
        %% update weights on new stimuli for initial calc of selectivity
        
        
        for cycle=1:p.numEncodingCycles
            
            [win_row, win_col, dist_mat] = ...
                findWinningNode(weights, input_mat, p.numInputDims(layer));
            
            
            p.winning(layer,grid,trial,1:2) = [win_row, win_col];
            
            %--------------------------------------------------------------
            %Calculate each unit's distance from winner and activation
            %--------------------------------------------------------------
            [f, selectivity, p, ~, ~] = ...
                calc_selectivity(win_row, win_col, dist_mat, p, p.numInputDims(layer));
            
            if cycle==1,  
                initial_selec(layer,grid) = selectivity;
            end
            
            % Update Weights
            weights = weights + f.*(input_mat-weights);  % update based on spire around winning node
            
            
        end 
        
        % load updated weights
        W(layer,:,:,1:p.numInputDims(layer),grid) = weights;
        
        
    end
end

end