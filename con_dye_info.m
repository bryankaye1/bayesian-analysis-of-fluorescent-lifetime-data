%%
% Here we are creating structures with the histogram (info.his) of long and short
% lifetimes control dyes (all the arriival times in one sample grouped into one
% histogram). Then we attach found/recorded parameters into the structre
% called info
%
clear
%%
loadname = 'll-his.mat';
dtemp =  load(loadname);
info = dtemp.data;
info.his = info.his;
info.int = 1.4*10^5/(8*10^7);
info.w2min = 3.5;
info.w2max = 4.5;
info.lifetime = 3.94;
save(loadname,'info');
%%

clear info tempdata
loadname = 'lm-his';
tempdata =  load(loadname);
info.his = tempdata.data;
info.int = 1.1*10^6/(8*10^7);
info.w2min = 3.5;
info.w2max = 4.5;
info.lifetime = 3.972;
save(loadname,'info');

clear info tempdata
loadname = 'lh-his';
tempdata =  load(loadname);
info.his = tempdata.data;
info.int = 4.3*10^6/(8*10^7);
info.w2min = 3.5;
info.w2max = 4.5;
info.lifetime = 4.062;
save(loadname,'info');


clear info tempdata
loadname = 'sl-his';
tempdata =  load(loadname);
info.his = tempdata.data;
info.int = 1.5*10^5/(8*10^7);
info.w2min = .1;
info.w2max = 1;
info.lifetime = .454;
save(loadname,'info');

clear info tempdata
loadname = 'sm-his';
tempdata =  load(loadname);
info.his = tempdata.data;
info.int = 1.5*10^6/(8*10^7);
info.w2min = .1;
info.w2max = 1;
info.lifetime = .447;
save(loadname,'info');

clear info tempdata
loadname = 'sh-his';
tempdata =  load(loadname);
info.his = tempdata.data;
info.int = 4.6*10^6/(8*10^7);
info.w2min = .1;
info.w2max = 1;
info.lifetime = .452;
save(loadname,'info');

