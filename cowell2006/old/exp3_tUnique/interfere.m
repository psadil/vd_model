function [ p, weights ] = interfere( p, weights )
%interfere simulates intereference during delay cycles

%%
% flag for calc_act_fast to indicate that we're using the learning rate for
% interference cycles
interefere = 1;


%%
for layer = 1:max(p.numLayers)
    for grid = 1:p.nGrids(layer),
        
        w = squeeze(weights(layer,:,:,1:p.numInputDims(layer),grid));
        
        for cycle=1:p.delayCycles,
                        
            %--------------------------------------------------------------            
            % generate interfering stim
            %--------------------------------------------------------------
            inp_mat = gen_limited_input(p.numInputDims(layer),p); 
            
            %--------------------------------------------------------------
            % find node that best matches generated stim
            %--------------------------------------------------------------
            [win_row, win_col, dist_mat] = findWinningNode(w, inp_mat, p.numInputDims(layer));
            
            
            %--------------------------------------------------------------
            % Calculate each unit's distance from winner, for use in
            % updating
            %--------------------------------------------------------------
            [f, ~] = calc_act_fast(win_row, win_col, dist_mat,layer,p, interefere);
             
            %--------------------------------------------------------------
            % Update Weights
            %--------------------------------------------------------------
            w = w + f.*(inp_mat-w);
            
            
        end
        %------------------------------------------------------------------
        % store new weights in original weight array
        %------------------------------------------------------------------
        weights(layer,:,:,1:p.numInputDims(layer),grid) = w;
    end
    

end 

end