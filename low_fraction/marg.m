% marg.m function
% This function takes in the likelihood and outputs the marginalized
% likilihoods (values and coordinates) along pr, w02, w2, and w1 dimensions.
% It fuction also outputs the coordinate (value and index) of the maximum 
% of the marginalized like
%
function[prest, prBesti, prBest, prestx, w1est, w1Besti,w1Best, w1estx, w2est,w2Besti,w2Best, w2estx,...
    w02est,w02Besti,w02Best, w02estx, extest,extBesti,extBest,extestx] =...
    marg(post,prstep,prmin,prmax,w1step,w1min,w1max,w2step,w2min,w2max,w02step,w02min,w02max,extstep,extmin,extmax)


prest = squeeze(sum(sum(sum(sum(post,1),2),3),5));  %Can replace with trapz's later.
prest = prest/sum(prest);
prestx = prmin:prstep:prmax; 
[~, prBesti]= max(prest);
prBest = (prBesti-1)*prstep+prmin; %Marginalize best estimate of pr


w1est = squeeze(sum(sum(sum(sum(post,2),3),4),5));
w1est = w1est/sum(w1est);
w1estx = w1min:w1step:w1max;
[~, w1Besti]= max(w1est);
w1Best = (w1Besti-1)*w1step+w1min; %Marginalize best estimate of w1


w2est = squeeze(sum(sum(sum(sum(post,1),3),4),5));
w2est = w2est/sum(w2est);
w2estx = w2min:w2step:w2max;
[~, w2Besti]= max(w2est);
w2Best = (w2Besti-1)*w2step+w2min; %Marginalize best estimate of w1


w02est = squeeze(sum(sum(sum(sum(post,1),2),4),5));
w02est = w02est/sum(w02est);
w02estx = w02min:w02step:w02max;
[~, w02Besti]= max(w02est);
w02Best = (w02Besti-1)*w02step+w02min; %Marginalize best estimate of w02


extest = squeeze(sum(sum(sum(sum(post,1),2),3),4));
extest = extest/sum(extest);
extestx = extmin:extstep:extmax;
[~, extBesti]= max(extest);
extBest = (extBesti-1)*extstep+extmin; %Marginalize best estimate of ext


end