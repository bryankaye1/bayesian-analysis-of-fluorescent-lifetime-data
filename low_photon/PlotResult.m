function PlotResult(method,saveornot,name,time,decay,nexpo,fitstart,fitend,pfit,pstd,pvec,margpost)

% PlotResult(method,saveornot,name,time,decay,nexpo,fitstart,fitend,pfit,pstd,pvec,margpost,post)
% This function make a figure showing FLIM curve, the best estimated model,
% and weighted residual
%
% method: 1 for LS, 3 for Bayes grid
% saveornot: boolean, save if 1, don't save if 0
% name: name of the file, if you want to save plot
% time: time axis of FLIM curve
% decay: y axis of FLIM curve
% nexpo: number of exponentials in FLIM model
% fitstart, fitend: indices defining ROI
% pfit: best estimated parameters
% pstd: errors (standard deviation) in best esimated parameters
% pvec: parameter grid (Bayes only)
% marpost: marginal posterior (Bayes only)

numtitle = 'off';
if saveornot == 0;
    name = 'Figure';
   numtitle = 'on';
end


if nexpo == 1
    Nparam = 3;
elseif nexpo == 2
    Nparam = 5;
end

totcounts = sum(decay(fitstart:fitend));

%weight on residual
nonzero_decay = decay;
nonzero_decay(decay==0)=1;
sigy = sqrt(nonzero_decay);
weight = 1./sigy(fitstart:fitend);


f1 = figure;
set(f1,'name',name,'numbertitle',numtitle)
subplot(4,1,[1;2;3])
semilogy(time,decay,'.r')
hold on;
y_hat = lm_decay_model(time,pfit,[nexpo,totcounts,fitstart,fitend]);
y_hat = y_hat(fitstart:fitend);
semilogy(time(fitstart:fitend),y_hat,'-k','linewidth',1.5);
yrange = ylim;
text(7,yrange(2)/2,['Counts : ',num2str(totcounts)]);

xlim([0,10])

y_dat = decay(fitstart:fitend);
weighted_residual = weight.*(y_dat-y_hat);    
subplot(4,1,4);
plot(time(fitstart:fitend),weighted_residual);
axis([0,10,-5,5])
hold on
line([0;10],[0;0],'Color','k');


if method == 3 %Bayes fit
    f2 = figure;
    set(f2,'name',name,'numbertitle',numtitle)
    for i = 1:Nparam
        subplot(Nparam,1,i);
        plot(pvec{i},margpost{i});
        hold on;
        str(1) = {['Avg : ',num2str(pfit(i))]};
        str(2) = {['Std : ',num2str(pstd(i))]};
        xrange = xlim;
        yrange = ylim;
        text(0.7*(xrange(2)-xrange(1))+xrange(1),0.8*yrange(2),str);
    end
end

if saveornot == 1
    saveas(f1,['SavedResult/',name,'_fit'],'fig');
    close(f1);
    if method ==3
        saveas(f2,['SavedResult/',name,'_margpost'],'fig');
        close(f2);
    end
end
