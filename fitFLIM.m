function [model] = fitFLIM(output,qv,ti,fignum)

T=12.58;
bins = output(qv).bins;
brem = output(qv).brem;
ga = output(qv).ga;
ga(1900:end) = 0; %Temporary line for messing around. REMOVE
binskeep = bins - brem;

s = T/bins:T/bins:T; %time vector used to generate PDF of signal exponential
wig = output(qv).wig; %wig = wig';
tfw =output(qv).tfw; %forward amount of time to remove
tbac =output(qv).tbac; %backwards amount of time to remove
p = output(qv).datahis;
w2 = output(qv).w2Best;
w02 = output(qv).w02Best;
w01 = output(qv).prBest;
r1l = output(qv).r1l;
r2l = output(qv).r2l;

[srem,erem,p2,wig2] = remove_bins(T,bins,tfw,tbac,p,wig); %returns data, wigs, and indeces of which bins to keep/remove

f2 = exp(-s/w2); %signal over one period
f2 = [f2 f2]; %signal over 2 consecutive periods
f2con = conv(f2,ga); %PDF after conv
f2bar = f2con(bins+1:2*bins); %pdf after mod-ing by laser period
f2h = f2bar(1:binskeep); %Keep only the appropriate bins
f2h = f2h/sum(f2h);

cf2h = fliplr(cumsum(fliplr(f2h))); %for int correction
f2h2 = f2h.*cf2h/sum(f2h.*cf2h); %fliplr(cumsum(fliplr(f2h))) is CDF.
f2hr = f2h*r1l + f2h2*r2l; %f2h and f2h2 are both normalized. Since r1l+r2l =1, f2hr is normalized

back = (1-w01-w02)/bins;
model = (f2hr + back).*wig;

model = model/sum(model);
pn = p/sum(p);
figure(fignum);clf; hold all; plot(pn, '.','MarkerSize',5, 'Color', 'b'); plot(model,'r'); title(ti)

end
