%param.m function
% This function finds the new bounds and grid point spacings 
% on the paramter search by analyzing the marginalized likelihood.
% param.m function
% It always adds sl/sr steps to the left/right of where the likelihood
% falls below the threshold, which is .01*max. It will never increase the
% bounds on the search space.
% error variable may not be doing anything.
function [w1minout, w1maxout, error] = param(w1est, thr, w1min, w1max, w1step, sl, sr)
sl = round(sl);
sr = round(sr);
error = 0;
if w1min-w1max~=0
    
    thr = thr*max(w1est);
    w1estm = w1est<thr;
    
    w1mini = find(w1estm==0, 1, 'first');
    w1maxi = find(w1estm==0, 1, 'last');
    
    w1minout = (w1mini-sl)*w1step+w1min;
    w1maxout = (w1maxi+sr)*w1step+w1min;
    
    if isempty(w1estm)
        w1minout = w1min;
        w1maxout = w1max;
    end
    
    %Make sure search space does not grow%
    if w1min > w1minout
        w1minout = w1min;
    end
    
    if w1max < w1maxout
        w1maxout = w1max;
    end
    
    %Im pretty sure these errors cannot be triggered, but just in case:
    if w1maxi >= length(w1estm)
        if w1max~=w1maxout
            error = pi;
        end
    end
    if w1mini <= 1
        if w1min~=w1minout 
        error = 1;
        end
    end
    
else   
    w1minout = w1min;
    w1maxout = w1max;  
end

    %Set errors if necessary%
    %     if w1estm(w1mini-1)> w1est(w1mini) || w1estm(w1maxi+1) > w1estm(w1maxi)
    %         error ='not monotonic about max';
    %     else
    %         error ='none';
    %     end


