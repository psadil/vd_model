function [inp_matrix] = gen_limited_input(numInputDims,p)
numRows = p.numRows;

nSimpleConj = (p.numInputDims_Caudal / p.nDimReps) ^ p.nStimFactors;
count=1;
avail_features = zeros(nSimpleConj, p.numInputDims_Caudal / p.nDimReps);
for inp1 = 1:p.numGrids_Caudal,
    for inp2 = 1:p.nStimFactors,
        avail_features(count,:) = [inp1 inp2];
        count=count+1;
    end
end


avail_features(avail_features==1)=0.05;
avail_features(avail_features==2)=0.35;
avail_features(avail_features==3)=0.65;
avail_features(avail_features==4)=0.95;


%%%%% Generate an input

if numInputDims == p.numInputDims_Caudal / p.nDimReps,
    segment = avail_features(randi(nSimpleConj),:);
    
    
    inp_matrix_temp = ones(numRows,numRows,p.numInputDims_Caudal / p.nDimReps);
    inp_matrix = zeros(size(inp_matrix_temp));
    for dim = 1:p.numInputDims_Caudal/p.nDimReps
        inp_matrix(:,:,dim) = segment(dim) * inp_matrix_temp(:,:,dim);
    end

    % Take the input x,y coords, and make each fill a slice
    % of a matrix the same size as m, so that they can be subtracted
elseif numInputDims == p.numInputDims_PRC / p.nDimReps
    for seg = 1:p.numGrids_Caudal ,
        segment(seg,:) = avail_features(randi(nSimpleConj),:);
    end
    
    inp_matrix_temp = ones(numRows,numRows,p.numInputDims_PRC / p.nDimReps);
    inp_matrix = zeros(size(inp_matrix_temp));
    for dim = 1:p.numInputDims_PRC / p.nDimReps
        inp_matrix(:,:,dim) = segment(dim) * inp_matrix_temp(:,:,dim);
    end
    
else
    fprintf('Invalid number of stimulus dimensions \n');
end

inp_matrix = reshape(repmat(reshape(inp_matrix,...
    [1,p.numRows,1,p.numRows,1,numInputDims]),[1,1,1,1,p.nDimReps]),...
    [p.numRows,p.numRows,numInputDims*p.nDimReps]);


end