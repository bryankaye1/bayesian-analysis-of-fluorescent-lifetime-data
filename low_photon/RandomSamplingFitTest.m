% This code generates FLIM curves composed of N number of photons by
% randomly sampling photons from the master FLIM data.
% 300 FLIM curves are sampled for each N, where N takes 11 different
% numbers from 50 to 10^5
% The sampled FLIM curves then analyzed using Bayesin FLIM approach to
% evaluate the accuracy and precision.
%
% Copyright (c) 2016, Tae Yeon Yoo
% all rights reserved
% tyoo@fas.harvard.edu

clc;
clear all;
close all;

addpath CONVNFFT_Folder

%% Choose the method you want to test:
method = 3;   % 3 for Bayes, and 1 for Least square fit
prior = 1;   % 1 for Constant prior


%% Load irf
loaded_irf = load('currentIRF.mat');
irf = loaded_irf.decay;
time_irf = loaded_irf.time;

dt_irf = time_irf(2) - time_irf(1);  %size of the time bin in irf
dt_macro = 25*10^-3;  %macrotime, in microsecond

%% Load master data
% fname = 'ex850nm_em542_27_035v_10x_int10sec_singlepoint_m1.asc';
% imported = dlmread(fname);

%% Read microtime and macrotime from the master data
% macroch = imported(:,1);
% microch = imported(:,2);
% macrot = (imported(:,1)-imported(1,1))*dt_macro;
% microt = imported(:,2)*dt_irf;
% flag = imported(:,4);

%% Set up parameters
t_start = 0.6;  
t_end = 9.4;  %t_start and t_end define the region of FLIM curve you wish to analyze
adc_ratio = 1; 
nexpo = 2;
Nparam = 5;

time = time_irf(adc_ratio:adc_ratio:length(time_irf));
% decay = hist(ceil(microch/adc_ratio),1:length(time));
% decay = decay';
% 
% totcounts = length(macroch);

dt = time(2)-time(1);  %size of time bin
fit_start = round(t_start/dt);    %define roi in terms of time bin #
fit_end = round(t_end/dt);


%% weight on residual
% nonzero_decay = decay;
% nonzero_decay(decay==0)=1;
% sigy = sqrt(nonzero_decay);
% weight = 1./sigy(fit_start:fit_end);



%% Random Sampling
Ncounts = [50;100;200;400;800;1600;3200;6400;12800;25600;100000];
Nsamples = 300;  %Nsamples samples per each Ncounts

sampled_decay = zeros(length(time),length(Ncounts),Nsamples);

if method == 3
    %for Bayes
    samplepavg = zeros(Nparam,length(Ncounts),Nsamples);
    samplepstd = zeros(Nparam,length(Ncounts),Nsamples);
    samplemargpost = cell(Nparam,length(Ncounts),Nsamples);
    samplepost = cell(length(Ncounts),Nsamples);
    samplemle = zeros(Nparam,length(Ncounts),Nsamples);
    
end

if method == 1
    %for LS
    samplepfit = zeros(Nparam,length(Ncounts),Nsamples);
    samplepstd = zeros(Nparam,length(Ncounts),Nsamples);
    sampleChisq = zeros(length(Ncounts),Nsamples);
    
end

% for i = 1:length(Ncounts)
%     for j = 1:Nsamples
%         sample = randsample(microch,Ncounts(i));
%         sampled_decay(:,i,j) = hist(ceil(sample/adc_ratio),1:length(time));
%     end
% end

%Save random sampling results. If random sampling is already done, import
%the result instead

%save('SavedResult/sampled_decay_counts50to100000_300samples.mat','sampled_decay');

imported_sample = load('sampled_decay_counts50to100000_300samples.mat');
imported_decay = imported_sample.sampled_decay;

for i = 1:length(Ncounts)
    for j = 1:Nsamples
        reshaped = reshape(imported_decay(:,i,j),adc_ratio,length(time));
        sampled_decay(:,i,j) = sum(reshaped,1)';
    end
end


%%

% long lifetime is fixed to 4.03, the ratio of short to long lifetime fixed
% to 0.12. (they are determined from FLIM curve with a lot of photons)
tau1assym = 4.03;
Eassym = 0.12;
Aassym = 0.997;

%% Bayes fit

prior = 1;   %constant prior
saveornot = 1;

% Set up the parameter grid
if nexpo == 1
    p_min = [-16,0.9938,0.1]';
    p_max = [-16,0.9938,8]';
    dp = [1,0.0025,0.01]';
elseif nexpo == 2
    p_min = [1,0.2,tau1assym,0.01,Eassym]';
    p_max = [1,1,tau1assym,1,Eassym]';
    dp = [1,0.0025,0.01,0.01,0.01]';
end

save('SavedResult/Bayessetting_fixedtau1_fixedE_constprior_300samples.mat','p_min','p_max','dp','prior');

% matlabpool;  %Do parallel computing to save time
h = waitbar(0,'Please wait');
for i = 1:length(Ncounts);
    for j = 1:Nsamples
        tic
        [pavg,sigp,pvec,post,margpost,mle] = bayes_fit(time,sampled_decay(:,i,j),dp,p_min,p_max,nexpo,prior,fit_start,fit_end);
        samplepavg(:,i,j) = pavg;
        samplepstd(:,i,j) = sigp;
        samplemargpost(:,i,j) = margpost;
        samplepost{i}{j} = post;
        samplemle(:,i,j) = mle;
        toc
        waitbar(((i-1)*Nsamples+j)/(Nsamples*length(Ncounts)),h);
%         PlotResult(method,saveornot,['Counts',num2str(Ncounts(i)),'_',num2str(j)],...
%             time,sampled_decay(:,i,j),nexpo,fit_start,fit_end,pavg,sigp,pvec,margpost,post)
    end
end

% matlabpool close;
save('SavedResult/BayesResult_fixedtau1_fixedE_constprior_300samples.mat','samplepavg','samplepstd','samplemargpost','samplemle','samplepost','pvec');
close(h);
    
