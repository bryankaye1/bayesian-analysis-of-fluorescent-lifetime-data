%This matlab script reads in the photon arrival times (units is bins) from 
%the BH .asc file and builds a histogram of arrival times (for lifetime fitting) 
%and builds a .mat file of the arrival bins.
clear
%% Convert Short-Lifetime Photon Master List
tic
%format is short lifetime - low ,med , high intensity. Then long lifetime -
%low, med, high intensity.
asc_filename = {'eryth-1p5e5_550-88nm-950ex-2000sec_m1.asc','eryth-1p5e6_550-88nm-950ex-200sec_m1.asc'...
    'eryth-4p6e6_550-88nm-950ex-100sec_m1.asc','coumarin_1p4e5_550-88nm_950ex-2000sec_m1.asc',...
    'coumarin_1p1e6_550-88nm_950ex_200sec_m1.asc','coumarin_4p3e6_550-88nm_950ex_200sec_m1.asc'};
mat_filename = {'sl','sm','sh','ll','lm','lh'};
for i = 1:length(asc_filename)
    %load in asc file
fileID = fopen(strcat('C:\Users\Bryan\Documents\MATLAB\data\2014-10-07\',asc_filename{i}));
data1 = textscan(fileID, '%f %f %f %f','HeaderLines', 19 );
fclose(fileID);

%convert to matlab vector and create a histogram
arrival_times_list = data1{2};
num_photons = length(arrival_times_list);%number of photons
edges = 1:4096; %edges for binning arrival times from asc file (asc file records which bin the photon arrived in, bins 1-4096)
arrival_time_histrogram = histc(arrival_times_list, edges);
fprintf('%s has %3.0f photons\n',mat_filename{i},num_photons);

%save
save(strcat(mat_filename{i},'-his'),'arrival_time_histrogram');
save(strcat(mat_filename{i},'-at'),'arrival_times_list');
end
toc