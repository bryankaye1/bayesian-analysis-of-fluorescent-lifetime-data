%% this code calculates the "actual" intensity of the dye given a measured intensity
%adjusting for intensity with these numbers (at least for high int) gave a near
% perfect changes in FRET fraction. see "arival time conversion notes" in
% evernote for list of intensities

lrate = 7.9491e7;%laser repetition rate given time window is 12.58ns %8*10^7;
deadtime = 10^-7; %deadtime of BH system
intensity.measured_rate = {1.5e5,1.5e6,4.6e6,1.4e5,1.1e6,4.3e6}; %intensity
intensity.identifier = {'sl','sm','sh','ll','lm','lh'};
% This converts measured intensity to actual intensity, taking deatime and
%two photons in one period into account
for i = 1:length(intensity.measured_rate)
    ranew = intensity.measured_rate{i}/(1-intensity.measured_rate{i}*deadtime); %adjusts rate to account for deadtime
    ranew2 = ranew/lrate; %nomarlized rate (in terms of photons per laser period)
    ranew3 = ranew2/(1-ranew2); %adjusts rate to account for only measuring 1 photon per period
    
    intensity.r1{i} = ranew3 / (ranew3 + ranew3^2); %calculate prob of 1 pho
    intensity.r2{i} = ranew3*ranew3 / (ranew3 + ranew3^2); %calculate prob of 2 pho
    fprintf('%s has an r1=%1.4f r2 is %1.4f\n',intensity.identifier{i}, intensity.r1{i},intensity.r2{i});
end
    save('intensity_mat','intensity')



% %% This part of the code calculates the conversion from measured intensity to
% %%actual intensity, but only accounting for perturbations from 2 photons
% %arriving in the same peiod (and only one being measured). It does not take
% %into account the effect of deadtime. This is the way we analyzed data the
% %first time we ran these experiments,
% 
% rm2 = measured_rate/lrate; %nomarlized rate (in terms of photons per laser period)
% raold = rm2/(1-rm2); %convert measured to actual intensity
% 
% ra1old = raold/(raold+raold^2); %calculate prob of 1 pho
% ra2old = (raold^2)/(raold+raold^2);%calculate prob of 2 pho
% %note denom is normalization constant. Since photons are independent, prob
% %of indepent events is just product of prob of events