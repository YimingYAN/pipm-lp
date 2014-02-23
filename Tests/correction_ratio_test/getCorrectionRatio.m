%% Function used to calculate the correction ratios
function [falsePrediction, missedPrediction, cr] = ...
    getCorrectionRatio(predictedActv,actualActv)

% actualActv is empty?
if isempty(actualActv)
    if ~isempty(predictedActv)
        falsePrediction = 1;
    else
        falsePrediction = 0;
    end
    
    missedPrediction = 0;
    cr = 1-falsePrediction;
else
    M = length(union(predictedActv,actualActv));
    
    falsePrediction =...
        length(setdiff(predictedActv,actualActv))/M;
    missedPrediction = ...
        length(setdiff(actualActv,predictedActv))/M;
    
    cr = 1-falsePrediction - missedPrediction;
end


end
