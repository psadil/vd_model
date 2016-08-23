function [p,weights] = pretrain(p)
%%%%%%%%%%%%%%%%%%%%% Perform pre-training of network %%%%%%%%%%%%%%%%

fprintf('\nVD_pretrain being executed...');
interfere = 0;
weights=zeros(p.numLayers,p.numRows,p.numRows,p.numInputDims(p.numLayers),max(p.nGrids));

trainGrid = ones(p.numLayers,max(p.nGrids));

%%
for layer = 1:p.numLayers
    
    nInpDims = p.numInputDims(layer);
    
    for grid = 1:p.nGrids(layer),
        fprintf('\n\nLayer number %d, Grid no. %d...\n(1==Caudal, 2==PRC)', layer, grid);
        
        if trainGrid(layer,grid)
            
            % initialize random weights for grids in this layer
            w = init_weights(p, layer);
            
            %--------------------------------------------------------------
            % begin training cycle of newly generated grid
            %--------------------------------------------------------------
            
            for batch = 1:p.nBatch
                
                % sum mse initialized at 0 for each batch
                smse = 0;
                for cycle=1:p.numTrainCycles,
                    
                    p.eta = ((p.numTrainCycles*(batch-1))+cycle)^(-p.A);		% Learning rate (how quickly weights are adapted)
                    p.G = p.sigma2 + (p.numRows/2 - p.sigma2)*...
                        ((p.numTrainCycles*(batch-1))+cycle)^(-p.B);		% Gaussian width parameter
                    
                    
                    %                     if cycle == 1
                    %                         fprintf('\nWithin pretrain, Cycle %d, G = %f, ETA = %f',...
                    %                             cycle, p.G, p.eta);
                    %                     end
                    
                    
                    % Generate input data that is p.nInpDims
                    % by p.nRows
                    [inp_mat] = gen_limited_input(nInpDims,p); %generate an input vector
                    
                    
                    %----------------------------------------------------------
                    % Find winning node
                    %----------------------------------------------------------
                    [win_row, win_col, dist_mat] = findWinningNode(w, inp_mat, nInpDims);
                    
                    
                    % Calculate each unit's distance from winner and resultant activation
                    [f, ~] = calc_act_fast(win_row, win_col, dist_mat,layer,p,interfere);
                    
                    
                    % Update Weights
                    w = w + f.*(inp_mat-w);
                    
                    smse = smse + min(dist_mat(:));
                end  % end of training cycles loop
                
                mse = smse / p.numTrainCycles;
                
                % test to see if mse has reached desired error
                if mse < p.mse;
                    trainGrid(layer,grid) = 0;
                    fprintf('\nWithin pretrain, batch %d, G = %f, ETA = %f, MSE = %d',...
                        batch, p.G, p.eta, mse);
                    break % training cycle loop
                end
                
                % display batch
                if mod(batch,100)==0,
                    fprintf('\nWithin pretrain, Cycle %d, G = %f, ETA = %f, MSE = %d',...
                        batch, p.G, p.eta, smse);
                end
                
                
            end % end of batch loop
            
            
            % load trained weights of this grid back into full weight array
            weights(layer,:,:,1:p.numInputDims(layer),grid) = w;
            
        else
            
            % if the grid has already reached an aceptable mse, don't train
            % it anymore
            continue
        end
        
    end % end of grid loop
    
end % end of layer loop

end