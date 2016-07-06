% Error check function:
% This function checks thebounds set on the parameter search space. 
%
%On loop 1 of post_int, this function checks to see if the marginalized 
% likelihood maximum is at one the edge of the search space. This does not
% apply to pr and w02 if the maximum is at 0 or 1, since we know we the
% actualy value for pr and w02 must lie between these bounds.
% produce an error if the bounds are already 
%
%On loop 2, it checks if the marginalized like is sufficiently far from
% bounds.  If the likelihood is larger than .1*max, then code will set:
% errX = 2 for lower bound error
% errX = 3 for upper bound error
% errX = 5 for upper and lower bound error (this can happen if marg is too
% wide). X stands for pr, w02, w1,or w2.

function [errPR, errw02, errw2, errw1] = errorcheck(prBesti,prest,prestx,prBest,...
                        w02Besti,w02est,w02estx,w02Best,...
                        w2Besti,w2est,w2estx,w1Besti,w1est,w1estx,l1)

errPR = 0;
errw02 = 0;
errw1 = 0;
errw2 = 0;

if l1==1
    if prBesti == length(prest) || prBesti == 1
        if prBest ~= 0 && prBest ~=1 && range(prest)~=0
            errPR = 1;
        end
    end    
    if w02Besti == length(w02est) ||w02Besti == 1
        if w02Best ~= 0 && w02Best ~=1 && range(w02est)~=0
            errw02 = 1;
        end
    end   
    if w2Besti == length(w2est) || w2Besti == 1
        if range(w2est)~=0
            errw2 = 1;
        end
    end   
    if w1Besti == length(w1est) || w1Besti == 1
        if range(w1est)~=0
            errw1 = 1;
        end
    end
end

if l1==2
    if range(prestx)~=0   
        if prestx(1)~=0 && (prest(1)/max(prest))>0.1
            errPR = 2+errPR;
        end
        if prestx(end)~=1 && (prest(end)/max(prest))>0.1
            errPR = 3+errPR;
        end
    end
           
    if range(w02estx)~=0   
        if w02estx(1)~=0 && (w02est(1)/max(w02est))>0.1
            errw02 = 2+errw02;
        end
        if w02estx(end)~=1 && (w02est(end)/max(w02est))>0.1
            errw02 = 3+errw02;
        end
    end
    
    if range(w2estx)~=0   
        if (w2est(1)/max(w2est))>0.1
            errw2 = 2+errw2;
        end
        if (w2est(end)/max(w2est))>0.1
            errw2 = 3+errw2;
        end
    end
    
    if range(w1estx)~=0   
        if (w1est(1)/max(w1est))>0.1
            errw1 = 2+errw1;
        end
        if (w1est(end)/max(w1est))>0.1
            errw1 = 3+errw1;
        end
    end
end

end