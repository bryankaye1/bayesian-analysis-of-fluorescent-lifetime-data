Readme for data presented in Figures 1 and 3 (low fraction regime). To recreate the figure 3, run "control_dye_plots_v2". To recreate the data in figure 1, run "fig1_plots.m".



To create data sets:
1) Run "con_dye_hist.m". This code creates 100 FLIM histograms for each intensity setting and save them as matfiles.




To analyze data: "data_prep..." packages all the necessary information needed to analyze data into one .mat file called matin. These create sequential matins, i.e. matin1, matin2... . Use "post_int_working.m" to read these matins and create matouts, .mat files that contain the matin information as well as posterior distributions.


0) Optional: Run "append_cdye_his.m". This code appends the intensity info and long lifetime search parameters to each dataset. Not necessary to reproduce data

2) Run "data_prep_cdyes_find_lifetimes.m". This code prepares data for finding the lifetimes of each dye taken every intensity.

3) Run "post_int_corr". This code analyzes the data and returns posterior distributions.

4) Run "read_control_dye_lifetimes" to save maximum a posteriori lifetimes to a matfile and print nice graphs.

5a) Run "data_prep_cdyes_find_phofraction.m" to prepare changes in photon fraction data
5b) Run "data_prep_cdyes_find_phofraction_changing_lifetimes.m" to prepare FLIM data witha fixed photon fraction but changing photon number.
 

6) Run "post_int_corr" to obtain posteriors on 5a and 5b. 

7) Run "control_dye_plots_v2" to make the plots of the results of the data in 5a & 5b (shown in Fig 3).

8a) To make the plots in 1, Run "data_prep_cdyes_50_50_phofraction.m" to prepare data for figure 1 posterior. 
8b) Run "post_int_corr" with the matin number from 8a.
8c) Run "fig1_plots.m"




Code not used in paper:
Intensity correction:
cdyes_int_calc
