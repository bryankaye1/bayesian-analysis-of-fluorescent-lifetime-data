%This function checks the loglike to make sure:
% (1) the entire space was searched ("errorsize" variable)
%       if there are NaNs in loglike, errorize = 111
%       if the whole search space wasn't searched, errorsize = 123
% 
%(2) there are non-zero probabilities in the likilihood ("errorinf")
%      if all probs are zero, then errorinf = 1
function [errorsize,errorinf] = checkloglike2(loglike)
errorsize = 0;
errorinf = 0;

if sum(sum(sum(sum(sum(isnan(loglike)))))) > 0
    disp('ERROR ERROR ERROR: NaNs in loglike');
    errorsize = 111;
end
b = loglike(loglike==-pi);
if b~=0
    disp('ERROR ERROR ERROR: Did not search whole space');
    errorsize =123;
end

if max(max(max(max(max(loglike)))))==-inf
    errorinf=1;
    disp('ERROR ERROR ERROR: All prob in like are 0');
end

end