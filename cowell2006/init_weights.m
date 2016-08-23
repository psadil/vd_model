function m_out = init_weights(p, layer)
% init_weights -- initializes weigths for a given grid in a given layer.
% Currently weights are initialized randomly. Other schemes exist if
% pretraining were important to simulation, see Kohonen's book.

% called by: pretrain
% calls: NA

% input:
%   p: experimental structure -- useful for input dims and num rows.
%   layer: which layer to pull 

% output:
%    m_out: matrix that goes out, filled with weights for a single grid in
%       a single layer

%%
m_out = rand(p.nRows,p.nRows,p.nInputDims(layer));


end


