function m_out = init_weights(p, layer)
% generate normalized weights

m_out_tmp = rand(p.numRows,p.numRows,p.numInputDims(layer));

m_out = m_out_tmp ./ ...
    repmat(squeeze(sqrt(sum((m_out_tmp).^2,3))),[1,1,numInputDims]);

end


