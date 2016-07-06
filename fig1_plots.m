
clear;
load('Y:\Users\bkaye\cluster\matout\matout7553.mat','-mat');
y = output(1,1,1).datahis;
figure(1);clf; plot(y);
set(gca, 'YScale', 'log','XTick',[0,2927/2,2927]);
axis([0 2927 50 2e4]);

%load('Y:\Users\bkaye\cluster\matout\matout7554.mat','-mat');
xmin = 0.4852; xmax = 0.486;
ymin = 0.49; ymax = 0.5;
post = output(1,1,1).posterior;
X = output(1,1,1).w02estx; %W02
Y = output(1,1,1).prestx; %W02
post = squeeze(post);
post = post / ( sum(sum(post)) * ((xmax-xmin)*(ymax-ymin))  );
figure(2); clf; contourf(X,Y,post);
set(gca,'XLim',[0.4852,0.486],'YLim',[0.49,0.5], 'XTick',.4852:.0002:.4860,...
  'YTick',.49:.002:.500); %xaxis is w02 and yaxis is w01
colorbar;

