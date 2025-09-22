function [ratioMat] = meanLogRatio(ROITable)

colNames = ROITable.Properties.VariableNames;
ratioMat = NaN(size(colNames,2), size(colNames,2), size(ROITable,1));
epsVal = eps;
for i_t = 1:size(ROITable)
    ratioTrial = NaN(size(colNames,2));
    for Ci = 1:numel(colNames)
        for Ri = 1:numel(colNames)
            p1 = ROITable.(colNames{Ri})(i_t);
            p2 = ROITable.(colNames{Ci})(i_t);

            % avoid log(0), avoid dividing by 0, happens when probability
            % is 0 or 1 so we clip it by the smalles number possible
            p1 = max(min(p1, 1 - eps), eps);
            p2 = max(min(p2, 1 - eps), eps);
            
            % Compute the log odds ratio (higher value - row element is higher than col element)
            logOddsRatio = log((p1 ./ (1 - p1)) ./ (p2 ./ (1 - p2)));
            
            % Store the mean log odds ratio in the matrix
            ratioTrial(Ri, Ci) = mean(logOddsRatio);
        end
    end
    ratioMat(:,:,i_t) = ratioTrial;
end
end