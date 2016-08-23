function [p] = VD_pretrain(p,rat)

global ROOT;

fprintf('\nVD_pretrain being executed...');

%%%%%%%%%%%%%%%%%%%%% Perform pre-training of network %%%%%%%%%%%%%%%%
for layer = 1:max(p.numLayers)
    
    nInpDims = p.numInputDims(layer);
    
    
    for grid = 1:p.nGrids(layer),
        fprintf('\n\nLayer number %d, Grid no. %d...\n(1==Caudal, 2==PRC)', layer, grid);
        
        
        % initialize random weights for grids in this layer
        w = VD_init_weights(p, layer);
        
        
        
        %------------------------------------------------------------------
        % begin training cycle of newly generated grid
        %------------------------------------------------------------------
        
        for cycle=1:p.numTrainCycles(layer),
            
            p.eta = cycle^(-p.A);		% Learning rate (how quickly weights are adapted)
            p.G = 0.5 + 10*cycle^(-p.B);		% Gaussian width parameter
            %%% Note: G<0.5 is boring because the Gaussian only covers one node
            
            if cycle == 1 || cycle == p.numTrainCycles(layer),
                fprintf('\nWithin pretrain, Cycle %d, G = %f, ETA = %f', cycle, p.G, p.eta);
            end
            
            %%% Generate input data that is p.nInpDims
            %%% by p.nRows ()
            [inp_mat] = VD_gen_limited_input(nInpDims,p); %generate an input vector
            
            %--------------------------------------------------------------
            % Find winning node
            %--------------------------------------------------------------
            [win_row, win_col, dist_mat] = findWinningNode(w, inp_mat, nInpDims);
            
            %------------------------------------------------------------------
            % Calculate each unit's distance from winner and resultant activation
            [f, ~] = VD_calc_act_fast(win_row, win_col, dist_mat,layer,p);
            

            %%% Update Weights
            w = w + f.*(inp_mat-w);
            
            
        end  % end of training cycles loop
        
        %------------------------------------------------------------------
        % Save this pretrained grid
        %------------------------------------------------------------------
        location = strcat(ROOT,'rats/rat', num2str(rat), '/pretrainedW__layer', num2str(layer), 'grid', num2str(grid),'.mat');
        save(location, 'w');
        
    end
        
end

end