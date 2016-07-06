%clear; %close all;
%This code reads the matouts provided by "post_int_corr" and creates mats
%with the maximum posteriori of the lifetimes for each dataset.

%Make sure to place the appropriate matin number in line 8
clear;
set(0,'DefaultTextInterpreter','none');
imin = 5033; imax = 5040; %Place the appropriate matin number here.
ind = 0;
lifetimes.identifier = {'sl','sm','sh','ll','lm','lh','sl-no_int_corr','ll-no_int_corr'};
for i = imin:imax
    ind = ind+1;
    load(strcat('Y:\Users\bkaye\cluster\matout\matout',num2str(i),'.mat'));
    ti = sprintf('%s w2:%1.3f',lifetimes.identifier{ind},w2Bestmat(1));
    fitFLIM(output,1,ti,2*ind-1);
    w2_stdev = sqrt(sum(output(1,1,1).w2est.*(output(1,1,1).w2estx.^2))...
    -sum(output(1,1,1).w2est.*output(1,1,1).w2estx).^2);
    figure(2*ind);clf; plot(output(1,1,1).w2estx,output(1,1,1).w2est);
    fprintf('lifetime for matin%1.0f is: %1.4f +/-%1.6f\n',i,w2Bestmat(1),w2_stdev);
    [ebinary,estring] = read_matout_errorcheck(output(1,1,1).error,output(1,1,1).errparam); %will print to command line if there is an error
    lifetimes.lifetime{ind} = w2Bestmat(1);
    lifetimes.stdev{ind} = w2_stdev;
    save('lifetimes.mat','lifetimes');
end
%%




