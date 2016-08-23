function [ p, weights ] = interfere( p, weights )
%delay_interfere simulates intereference during delay cycles
%   Detailed explanation goes here

fprintf('\n%d interference cycles being executed...', p.delayCycles(1));
interefere = 1;

for layer = 1:max(p.numLayers)
    for grid = 1:p.nGrids(layer),
        
        % put the variable 'weights' into the format previously accepted by the model.
        w = squeeze(weights(layer,:,:,1:p.numInputDims(layer),grid));
        
        for cycle=1:p.delayCycles(1),
                        
                        
            inp_mat = gen_limited_input(p.numInputDims(layer)/p.nDimReps,p); %generate an input vector
            
            %--------------------------------------------------------------
            % Find winning node
            %--------------------------------------------------------------
            [win_row, win_col, dist_mat] = findWinningNode(w, inp_mat, p.numInputDims(layer));
            
            
            %--------------------------------------------------------------
            % Calculate each unit's distance from winner, for use in
            % updating
            [f, ~] = calc_act_fast(win_row, win_col, dist_mat,layer,p, interefere);
                        
            %%% Update Weights
            w = w + f.*(inp_mat-w);
            
            
        end
        weights(layer,:,:,1:p.numInputDims(layer),grid) = w;
    end
    

end 

end

