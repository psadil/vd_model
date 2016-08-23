function [inp_matrix] = gen_randInp(p, nInputDims)
% genRandInp -- generates a sample stimulus for pretraining phase.
% All abstract features are fair game as stimulus. Needs to pull from
% entire stim space, since pretraining should allow for both global and
% local order.

% called by: pretrain
% calls: permn

% input:
%   numInputDims: dims attended to by current grid
%   p: experimental structure

% output
%   inp_matrix: sample stimulus presented to grid

%%

% simple conjunctions available to a given network. features are 1-d,
% simple conjunctions are defined as what gets fed to caudal layer. So,
% since caudal layer currently sees 2-d inputs, simple conjunctions are
% composed of 2 features. There are 4 possible levels of each feature, so
% there are 4^2 possible simple conjunctions. This will define number of
% rows of the matrix, avail_features

avail_features = permn(p.features, p.nInputDims_Caudal);

% generate every necessary simple feature conjunction segment
% how many segments needed depends on to which layer is receiving the input
% each seg defines a simple conjunction (and so has a length equal to
% number of dims input to caudal layer). The rows of 'segment' is equal to
% the number of simple conjunctions looked at by the grid. For caudal, this
% is just 1. For PRC, this is 4.
for seg = 1:nInputDims / p.nInputDims_Caudal
    segment(seg,:) = avail_features(randi(size(avail_features,1)),:);
end

% create stim out of those segments. That is, node must have a vector
% associated with it, where the vector contains all simple conjunctions
% defined in segment.
inp_matrix_temp = ones(p.nRows,p.nRows,nInputDims);
inp_matrix = zeros(size(inp_matrix_temp));
for dim = 1:nInputDims
    inp_matrix(:,:,dim) = segment(dim) * inp_matrix_temp(:,:,dim);
end


end