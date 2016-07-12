function [srem,erem,p2,wig2] = remove_bins(T,bins,tfw,tbac,p,wig)

%"bins" is the number of bins that make up 12.58ns
%%Changes vectors into the correct shape for removing time bins
if isrow(p)==0
    p =p';
end
if isrow(wig)==0
    wig =wig';
end

%%Find max
[~, im] = max(p);

tpb = T/bins; %Time interval per bin
bfw = round(tfw/tpb); %Number of bins to remove after the peak
bbac = round(tbac/tpb); %Number of bins to remove before the peak

%%Remove bins around max
srem = im - bbac; %index to start removing bins
erem = im + bfw; %index to stop removing bins

p2 = [p(1:srem-1),p(erem:end)]; %data after removing bins
wig2 = [wig(1:srem-1),wig(erem:end)]; %wiggles after removing bins


end

