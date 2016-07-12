
clear;
load('C:\Users\Bryan\Documents\MATLAB\FLIM\cdyes\cdata_low_100rep.mat','-mat');
y = data(1,1,1).his;
figure(1); clf;
plot(y);%(928:3855));
%set(gca, 'YScale', 'log');

new_vec = [];
istart = 1;
for i = 1:length(y)
    new_vec(istart:istart+y(i)-1) = i;
istart = istart + y(i);
end

edges = 0:length(y);
n_samples = round(length(new_vec)*.01);
indlong = randi(length(new_vec),1,n_samples); %generate index of photons
samples = new_vec(indlong); % pull photons
sample_his = histc(samples, edges); %histogram data
figure(2);clf; plot(sample_his);% set(gca, 'YScale', 'log');