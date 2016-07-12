% This code plots the results of random sampling test, genarated by the
% code "RandomSamplingFitTest.m". If you haven't run this code, do this
% first.


addpath CONVNFFT_Folder\

close all;
clear all;
clc;

Nparam = 5;

Ncounts = [50;100;200;400;800;1600;3200;6400;12800;25600;100000];
Nsamples = 300;  %Nsamples samples per each Ncounts

%fontsize
fs = 15;

fitstart = 246;
fitend = 3848;

%% Import Result

import_fixtau1fixE = load('SavedResult/BayesResult_fixedtau1_fixedE_constprior_300samples.mat');

samplepavg_fixtau1fixE = import_fixtau1fixE.samplepavg;
samplepstd_fixtau1fixE = import_fixtau1fixE.samplepstd;
samplemle_fixtau1fixE = import_fixtau1fixE.samplemle;
samplepost_fixtau1fixE = import_fixtau1fixE.samplepost;
samplepvec = import_fixtau1fixE.pvec;

%% convert A,f to wf, wnf, wb

wf = zeros(11,300);
wnf = zeros(11,300);
time = (1:4096)'*10/4096;
for i = 1:11
    for j = 1:300
        pvec = samplepavg_fixtau1fixE(:,i,j);
        
        fdecaymodel = lm_decay_model(time,[pvec(1:2);pvec(3)*pvec(5)],[1,-1,fitstart,fitend]);
        nfdecaymodel = lm_decay_model(time,[pvec(1:2);pvec(3)],[1,-1,fitstart,fitend]);
        decaymodel = lm_decay_model(time,pvec,[2,-1,fitstart,fitend]);
        
        fnum = sum(fdecaymodel(fitstart:fitend))*(1-pvec(4));
        nfnum = sum(nfdecaymodel(fitstart:fitend))*pvec(4);
        denom = sum(decaymodel(fitstart:fitend));
        
        wf(i,j) = fnum/denom;
        wnf(i,j) = nfnum/denom;
    end
end
wb = 1-wf-wnf;

%%

wf_mle = zeros(11,300);
wnf_mle = zeros(11,300);
time = (1:4096)'*10/4096;
for i = 1:11
    for j = 1:300
        pvec = samplemle_fixtau1fixE(:,i,j);
        
        fdecaymodel = lm_decay_model(time,[pvec(1:2);pvec(3)*pvec(5)],[1,-1,fitstart,fitend]);
        nfdecaymodel = lm_decay_model(time,[pvec(1:2);pvec(3)],[1,-1,fitstart,fitend]);
        decaymodel = lm_decay_model(time,pvec,[2,-1,fitstart,fitend]);
        
        fnum = sum(fdecaymodel(fitstart:fitend))*(1-pvec(4));
        nfnum = sum(nfdecaymodel(fitstart:fitend))*pvec(4);
        denom = sum(decaymodel(fitstart:fitend));
        
        wf_mle(i,j) = fnum/denom;
        wnf_mle(i,j) = nfnum/denom;
    end
end
wb_mle = 1-wf_mle-wnf_mle;

%%
truef = 0.67;
truepop = 0.997;
shift = 1;
tau1 = 4.03;
E = 0.12;

truepvec = [shift;truepop;tau1;truef;E];
fdecaymodel = lm_decay_model(time,[truepvec(1:2);truepvec(3)*truepvec(5)],[1,-1,fitstart,fitend]);
nfdecaymodel = lm_decay_model(time,[truepvec(1:2);truepvec(3)],[1,-1,fitstart,fitend]);
decaymodel = lm_decay_model(time,truepvec,[2,-1,fitstart,fitend]);

fnum = sum(fdecaymodel(fitstart:fitend))*(1-truepvec(4));
nfnum = sum(nfdecaymodel(fitstart:fitend))*truepvec(4);
denom = sum(decaymodel(fitstart:fitend));

truewf = fnum/denom;
truewnf = nfnum/denom;


%% Calculates the lifetime error, SEM, and average estimated standard deviation
wferr_fixtau1fixE = mean(wf,2)-truewf;
semwf_fixtau1fixE = std(wf,0,2)/sqrt(Nsamples-1);

wfmleerr_fixtau1fixE = mean(wf_mle,2)-truewf;
stderrwfmle_fixtau1fixE = std(wf_mle,0,2)/sqrt(Nsamples-1);


%%
close all
hfig_err = figure
errorbar(Ncounts,mean(wf,2),semwf_fixtau1fixE,'s','MarkerSize',10);

hold on
errorbar(Ncounts,mean(wf_mle,2),stderrwfmle_fixtau1fixE,'.','MarkerSize',30);

xlim([10,10^6])
ylim([0.05,0.2])
set(gca,'XScale','log')
legend('Posterior Mean','Posterior Mode')
line([10,10^6],[0,0],'Color',[0.5,0.5,0.5],'LineWidth',1)
xlabel('Photon Counts','FontSize',fs)
ylabel('$\hat{w}_f$','FontSize',15,'interpreter','latex')
set(gca,'FontSize',fs)
box off

print(hfig_err,'FractionErr_300samples_BayesOnly_Photon','-depsc2')
%%
close all
hfig_std = figure;
hstd = scatter(Ncounts,semwf_fixtau1fixE,100,'o','filled');
hold on

fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,-inf],...
               'Upper',[Inf,Inf],...
               'StartPoint',[0.05 -0.7],'Display','iter',...
               'TolFun',1E-11,'TolX',1E-11);
ft = fittype('a*x^b','options',fo);
ff = fit(Ncounts(:),semwf_fixtau1fixE,ft,'Exclude',1:4);
ci = confint(ff);

xlim([30,2*10^5])
xl = xlim;
xx = linspace(xl(1),xl(2),100);
plot(xx,feval(ff,xx),'Color',[0.5,0.5,0.5],'LineWidth',2)

set(gca,'XScale','log','YScale','log')
xlabel('Photon Counts','FontSize',fs)
ylabel('$Std(\hat{w}_f)$','FontSize',fs,'interpreter','latex')
set(gca,'FontSize',fs,'XTick',[1e2,1e3,1e4,1e5])


print(hfig_std,'FractionSampleStd_300samples_BayesOnly_photon','-depsc2')




