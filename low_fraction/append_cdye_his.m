%%
clear
%%

clear data
loadname = 'll-his.mat';
dtemp =  load(loadname);
data = dtemp.data;
data.his = data.his;
data.int = 1.4*10^5/(8*10^7);
data.w2min = 3.5;
data.w2max = 4.5;
save(loadname,'data');
%%


clear data
loadname = 'lm-his';
tempdata =  load(loadname);
data.his = tempdata.hislong;
data.int = 1.1*10^6/(8*10^7);
data.w2min = 3.5;
data.w2max = 4.5;
save(loadname,'data');

clear data
loadname = 'lh-his';
tempdata =  load(loadname);
data.his = tempdata.hislong;
data.int = 4.3*10^6/(8*10^7);
data.w2min = 3.5;
data.w2max = 4.5;
save(loadname,'data');


clear data
loadname = 'sl-his';
tempdata =  load(loadname);
data.his = tempdata.hisshort;
data.int = 1.5*10^5/(8*10^7);
data.w2min = .1;
data.w2max = 1;
save(loadname,'data');

clear data
loadname = 'sm-his';
tempdata =  load(loadname);
data.his = tempdata.hisshort;
data.int = 1.5*10^6/(8*10^7);
data.w2min = .1;
data.w2max = 1;
save(loadname,'data');

clear data
loadname = 'sh-his';
tempdata =  load(loadname);
data.his = tempdata.hisshort;
data.int = 4.6*10^6/(8*10^7);
data.w2min = .1;
data.w2max = 1;
save(loadname,'data');