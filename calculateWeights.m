function weights = calculateWeights(labels_cell)
    weights = ones(size(labels_cell));
    countc = 0;
    countp = 0;
    
    for i = 1:numel(labels_cell)
        if labels_cell{i} == 'controls'
            countc = countc + 1;
        elseif labels_cell{i} == 'patients'
            countp = countp + 1;
        end
    end

    for l = 1:numel(labels_cell)
        if labels_cell{l} == 'patients'
            weights(l) = countc/countp;
        end
    end 
end