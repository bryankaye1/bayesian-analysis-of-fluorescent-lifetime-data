%function [hislong,hisshort, data,dshortt] = read_spc_fn_v2(filename,mexp)

%% Get Long-Lifetime Photon Master List
%Long lifetime dye and short lifetime dyes were seperately measured at
%different intensities and compiled into histograms.
%
% type =1 for changing photon number, fixing photon fraction
% type =0 for changing fret fraction, fixing photon number to 5e7
%The first letter refers to the lifetime of the dye:
% "l" is for long
% "s" is for short
%The next letter refers to the intensity at which the data was taken
% "l" is for low
% "m" is for medium
% "h" is for high
% "-at" stands for "arrival times". Here we recorded actual arrival times
% (before histogram-ing). We draw photons from this distribution. 
% 
% "-his" stands for histogram. After drawing photons from the "-at"
% short and long lifetime distributions, photons are binned into histograms
% and filenames saved with the "-his" suffix.
%
% e.g. "lh-at" is the long-lifetime, high-intensity arrival time data
%      "sl-at" is the short-lifetime, low-intensity arrival time data      
% 
% Note to self: clean up code by using cells.
%
%
%% Combine Photons to make distributions
clear
type=1;

    file_names{1}{1} = 'll-at';
    file_names{1}{2} = 'sl-at';
    file_names{2}{1} = 'lm-at';
    file_names{2}{2} = 'sm-at';
    file_names{3}{1} = 'lh-at';
    file_names{3}{2} = 'sh-at'; 
    edges = 1:4096;
    
if type ==1  
    for oloop =1:3
        clear tempdlong.dlong tempdshort.dshort data dlong dshort;
        %load in appropriate arrival time data
        tempdlong = load(file_names{oloop}{1});
        dlong = tempdlong.arrival_times_list;
        tempdshort = load(file_names{oloop}{2});
        dshort = tempdshort.arrival_times_list;
        
        %Build histograms
        yimax = 8;
        repeatmax = 100;       
        data(yimax,repeatmax).his=pi;
        data(yimax,repeatmax).nplong = pi;
        data(yimax,repeatmax).npshort = pi; 
        data(yimax,repeatmax).hislong = pi;
        data(yimax,repeatmax).hisshort = pi;     
        tic
        for yi = 1:yimax
            fprintf('yi is %f\n',yi);
            y = 2^(-yi)*10^8;
            for repeat = 1:repeatmax  
                photon_fraction = 2^(-7); %This sets the fraction of short lifetime photons to toal photons
                
                nplong = round(y*(1-photon_fraction)); %number of long photons
                indlong = randi(length(dlong),1,nplong); %generate index of photons
                long = dlong(indlong); % pull photons
                hislong = histc(long, edges); %histogram data
                
                npshort = round(y*photon_fraction); % Number of short photons
                indshort = randi(length(dshort),1,npshort); %index for random sampling from short dist
                short = dshort(indshort); %pull photons
                hisshort = histc(short, edges); %histogram data
                
                data(yi,repeat).his=hisshort+hislong; % 'Mix' photons together and save as histogram analyze
                data(yi,repeat).nplong = nplong; %Save number of long lifetime photons in the matfile
                data(yi,repeat).npshort = npshort; %Sav number of short lifetime photons in the matfile                
                data(yi,repeat).hislong = hislong; %Save long lifetime histogram
                data(yi,repeat).hisshort = hisshort; %Save short lifetime histogram
                
            end
        end
        toc

        if oloop ==1
            save('cdye_ratio128x_lowInt','data'); %Name the file depending on the fraction
        elseif oloop ==2
            save('cdye_ratio128x_medInt','data');
        elseif oloop ==3
            save('cdye_ratio128x_hiInt','data');
        end    
    end
end

if type ==0    
    for oloop =1:3
        clear tempdlong.dlong tempdshort.dshort data dlong dshort;
        %Load in appropriate intensity data
        tempdlong = load(file_names{oloop}{1});
        dlong = tempdlong.arrival_times_list;
        tempdshort = load(file_names{oloop}{2});
        dshort = tempdshort.arrival_times_list;
        
        %Build histograms
        y = 5*10^7; %Set total number of photons   
        xmax = 18; %Sets the number of different photon fractions (from 1/2 to 1/2^xmax)
        repeatmax = 100; %Number of data sets.
        data(xmax,repeatmax).his=pi;
        data(xmax,repeatmax).nplong = pi;
        data(xmax,repeatmax).npshort = pi;        
        data(xmax,repeatmax).hislong = pi;
        data(xmax,repeatmax).hisshort = pi;
        
        tic
        for x = 1:xmax
            fprintf('x is %f\n',x);
            for repeat = 1:repeatmax
                
                photon_fraction = 2^(-x); %change photon fraction from 1/2 to 1/2^18               
                nplong = round(y*(1-photon_fraction)); %number of long photons
                indlong = randi(length(dlong),1,nplong); %Randomly sample index of photons
                long = dlong(indlong);%Randomly sampled long lifetime photons
                hislong = histc(long, edges); %histogram of long lifetime photons
                
                npshort = round(y*photon_fraction); % Number of short photons
                indshort = randi(length(dshort),1,npshort); %index for random sampling from short dist
                short = dshort(indshort); %Randomly pulled short lifetime photons
                hisshort = histc(short, edges); %histogram of short lifetime photons
                
                data(x,repeat).his=hisshort+hislong; %Add photon histograms
                data(x,repeat).nplong = nplong; %number of long lifetime phtons
                data(x,repeat).npshort = npshort; %number of short lifetime photons
                
                data(x,repeat).hislong = hislong;% histogram of long lifetime photons
                data(x,repeat).hisshort = hisshort;% histogram of short lifetime photons
                
            end
        end
        toc
        
        if oloop ==1
            save('cdata_low_100rep', 'data');
        elseif oloop ==2
            save('cdata_med_100rep', 'data');
        elseif oloop ==3
            save('cdata_hi_100rep', 'data');
        end
        
    end
    toc
end

