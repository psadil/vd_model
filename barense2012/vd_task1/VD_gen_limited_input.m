function [inp_matrix] = VD_gen_limited_input(numInputDims,p)
numRows = p.numRows;


avail_features = permn(p.features,nDimsCaudal);



if numInputDims == p.numInputDims_Caudal,
    segment = avail_features(randi(p.components),:);
    
    
    inp_matrix_temp = ones(numRows,numRows,p.numInputDims_Caudal);
    inp_matrix = zeros(size(inp_matrix_temp));
    for dim = 1:p.numInputDims_Caudal
        inp_matrix(:,:,dim) = segment(dim) * inp_matrix_temp(:,:,dim);
    end
    
    
    
elseif numInputDims == p.numInputDims_PRC
    for seg = 1:p.numGrids_Caudal ,
        segment(seg,:) = avail_features(randi(p.components),:);
    end
    
    inp_matrix_temp = ones(numRows,numRows,p.numInputDims_PRC);
    inp_matrix = zeros(size(inp_matrix_temp));
    for dim = 1:p.numInputDims_PRC
        inp_matrix(:,:,dim) = segment(dim) * inp_matrix_temp(:,:,dim);
    end
    
else
    fprintf('Invalid number of stimulus dimensions \n');
end

end