function [inp_matrix] = VD_gen_limited_input(numInputDims,p)

% gather list of all possible combinations of simple feature conjunctions
% current stim structure: should be 16, 2D simple conjunctions
nSimpleConj = p.numInputDims_Caudal ^ p.nStimFactors;
count=1;
avail_features = zeros(nSimpleConj, p.numInputDims_Caudal);
for inp1 = 1:p.numGrids_Caudal,
    for inp2 = 1:p.nStimFactors,
        avail_features(count,:) = [inp1 inp2];
        count=count+1;
    end
end

% assign values to each of the levels of all the input dims
avail_features(avail_features==1)=0.05;
avail_features(avail_features==2)=0.35;
avail_features(avail_features==3)=0.65;
avail_features(avail_features==4)=0.95;



% generate every necessary simple feature conjunction segment
% how many segments needed depends on to which layer is receiving the input
for seg = 1:numInputDims / p.numInputDims_Caudal
    segment(seg,:) = avail_features(randi(size(avail_features,1)),:);
end

% create stim out of those segments
inp_matrix_temp = ones(p.numRows,p.numRows,numInputDims);
inp_matrix = zeros(size(inp_matrix_temp));
for dim = 1:numInputDims
    inp_matrix(:,:,dim) = segment(dim) * inp_matrix_temp(:,:,dim);
end



end