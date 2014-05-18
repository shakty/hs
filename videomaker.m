%% Perform Cluster Analysis on the dump
tic
close all
clc
clear

path(path,'util/');

% Change default axes fonts.
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 14)

% Change default text fonts.
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 14)


% R = 0;

DUMPDIR = '/mnt/tmp/dump/EXP_BUG/';

MYDIR = 'attrLinear_nav_rndseeds_rndseq_tm_R0_n100_fv0';
%MYDIR = 'attrLinear_nav_rndseeds_rndseq_tm_Alpha1_n100_fv0';
MYDIR = 'attrZero_nav_rndseeds_rndseq_tm_R0_n100_fv0';
MYDIR = 'attrHard_nav_rndseeds_rndseq_tm_R0_n100_fv0';
MYDIR = 'attrMillean_nav_rndseeds_rndseq_tm_R0_n100_fv0';

% SIZE 10

DUMPDIR = '/mnt/tmp/dump/SIZE10/';

% R CLEAN
MYDIR = 'attrZero_nav_rndseeds_rndseq_tm_RClean_n100_fv0_s10';
MYDIR = 'attrLinear_nav_rndseeds_rndseq_tm_RClean_n100_fv0_s10';
MYDIR = 'attrMillean_nav_rndseeds_rndseq_tm_RClean_n100_fv0_s10';
% MYDIR = 'attrHard_nav_rndseeds_rndseq_tm_RClean_n100_fv0_s10'; (TODO);


% EXP_BUG

DUMPDIR = '/mnt/tmp/dump/EXP_BUG/';

% R CLEAN
MYDIR = 'attrLinear_nav_rndseeds_rndseq_tm_Rclean_n100_fv0';
MYDIR = 'attrMillean_nav_rndseeds_rndseq_tm_RClean_n100_fv0';
MYDIR = 'attrHard_nav_rndseeds_rndseq_tm_RClean_n100_fv0';
MYDIR = 'attrZero_nav_rndseeds_rndseq_tm_RleftClean_n100_fv0';


% EXP_BUG

DUMPDIR = '/mnt/tmp/dump/NAV/';

% R CLEAN
MYDIR = 'attrLinear_nav_rndseeds_rndseq_tm_Rclean_n100_fv0_s1';




% SIGMA
SIGMA = '_s1';


% DUMP DIR
dumpDir = [DUMPDIR MYDIR '/' MYDIR SIGMA '/'];
% FILE
myFile = '1-1.mat';

videoSubDir = [ MYDIR '/' ];
VIDEODIR = '/home/stefano/hs/videos/LAST/';
SAVE_VIDEO = 0;
SHOW_POTENTIAL = 0;
plottype = 3;
SINGLE = 1;

% PLOT TYPE
% plot_cross = 0;
% plot_number = 1;
% plot_number_color = 2;
% plot_arrow = 3;

% SINGLE
DUMPDIR = '/mnt/tmp/dump/NAV/';
simName = 'attrLinear_nav_rndseeds_rndseq_tm_RClean_n100_fv0_s1/attrLinear_nav_rndseeds_rndseq_tm_RClean_n100_fv0_s1_s2/';

DUMPDIR = '/home/stefano/hs/test/'; 
simName = 'VAGAIN-2013-12-18-15-11/'; % noise angular;
%simName = 'VAGAIN-2013-12-18-15-35/'; % noise on p;
%simName = 'VAGAIN-2013-12-18-15-56/'; % noise on p (R0.01 ALPHA0.1)
%simName = 'VAGAIN-2013-12-18-16-8/'; % no social influence noise on p
%simName = 'VAGAIN-2013-12-18-16-26/'; % no social influence noise on a small v, small noise (0.01)
%simName = 'VAGAIN-2013-12-18-16-32/'; % no social influence noise on a high v, small noise (0.01)
%simName = 'VAGAIN-2013-12-19-10-27/'; % no social influence noise on v reduced bouncing (on actual move)
simName = 'VAGAIN-2013-12-19-11-2-a/'; % social, noise a + noise on p (sigma x10)
%simName = 'VAGAIN-2013-12-19-11-10/'; % social, noise a + noise on p (sigma both 0.1)
% simName = 'VAGAIN-2013-12-19-11-14/'; % social, noise a (sigma 0.1) -> agents forms some clusters then go the truth
% simName = 'VAGAIN-2013-12-19-11-21/'; % social, noise a (sigma 0.1) + repulsion 1.5 -> agents forms a circle very large around
% simName = 'VAGAIN-2013-12-19-11-34/'; % social, noise a (sigma 0.1) (like above reduce velocity / 10) + repulsion 1.5
simName = 'VAGAIN-2013-12-19-11-37/'; % social, noise a (sigma 0.1) (like first reduce velocity / 10) B = 1.2 + repulsion 1.5 -> circle has a small radius
%simName = 'VAGAIN-2013-12-19-11-46/'; % social, noise a (sigma 0.01) -> agents go quickly in the middle three cluster forms very compact, very close to truth
%simName = 'VAGAIN-2013-12-19-11-54/'; % social, noise a (sigma 0.01) + noise also on average_velocity. -> agents increase / decrease velocity. three clusters forms going around. cool
%simName = 'VAGAIN-2013-12-19-11-58/'; % social, noise a (sigma 0.01) + noise also on average_velocity. -> agents increase / decrease velocity. clusters forms going around. cool
%simName = 'VAGAIN-2013-12-19-12-6/'; % like above + noise on p (like before a bit more spread (anyway noise on p is too small))
% simName = 'VAGAIN-2013-12-19-12-11/'; % like above + higher v (x10) -> more clusters and more spread

%simName = 'VAGAIN-2013-12-19-12-17/'; % like above + sigma on p (x10) -> many clusters less cohesive
% simName = 'VAGAIN-2013-12-19-12-22/'; % like above -> R much larger = .1 -> a big cluster and a small one
% simName = 'VAGAIN-2013-12-19-12-27/'; % no truth, R= 0.03, ALPHA = 0.1 velocity = 1 (like above), tau 1. some clustering. need to check with the default impl.
% simName = 'VAGAIN-2013-12-19-12-33/'; % Millean, R= 0.03, ALPHA = 0.1 velocity = 1 (like above), tau 1. some clusters in the radius, similar to standard settings
simName = 'VAGAIN-2013-12-20-11-13/'; % Linear nav np
% simName = 'VAGAIN-2013-12-20-11-23/'; % Small repulsion B = 1.001 navnp

simName = 'NAVPEPSILON-2013-12-21-12-51/'; % Epsilon 0.1 Sigma 0.01
simName = 'NAVPEPSILON-2013-12-21-13-3/'; % Epsilon 0.1 Sigma 0.05 -> clusters go towards the middle and merge
simName = 'NAVPEPSILON-2013-12-21-13-5/'; % Epsilon 0.5 Sigma 0.05 -> too noisy
simName = 'NAVPEPSILON-2013-12-21-13-10/'; % Epsilon 0.4 Sigma 0.04 -> too noisy
simName = 'NAVPEPSILON-2013-12-21-13-13/'; % Epsilon 0.3 Sigma 0.03 -> kind of too noisy
simName = 'NAVPEPSILON-2013-12-21-13-15/'; % Epsilon 0.2 Sigma 0.02 -> kind of OK for clusters
simName = 'NAVPEPSILON-2013-12-21-13-19/'; % Epsilon 0.1 Sigma 0.01 -> very good for clusters
simName = 'NAVPEPSILON-2013-12-21-13-23/'; % Epsilon 0.1 Sigma 0 -> very good for clusters (stable)
%simName = 'NAVPEPSILON-2013-12-21-13-26/'; % Epsilon 0 Sigma 0 -> perfect stable clusters
%simName = 'NAVPEPSILON-2013-12-21-13-30/'; % Epsilon 0.1 Sigma 0.01 R=0.1 -> 1 cluster near to truth
%simName = 'NAVPEPSILON-2013-12-21-13-32/'; % Epsilon 0.01 Sigma 0.01 -> 1 clusters very cohesive

simName = 'VELTEST-2013-12-28-12-39/'; % test with velocity: initV = 0.5 -> still some clustering both feeble 
% simName = 'VELTEST-2013-12-28-12-43/'; % like above, initV = 0.2 -> 1 big cluster and a few outsiders
% simName = 'VELTEST-2013-12-28-12-45/'; % like above, initV = 0.1 -> 1 big cluster
% simName = 'VELTEST-2013-12-28-12-47/'; % like above, initV = 1.5 -> many small clusters, some bouncing
simName = 'VELTEST-2013-12-28-12-49/'; % like above, initV = 2 -> more bouncing, many more spread out clusters
% simName = 'VELTEST-2013-12-28-12-50/'; % like above, initV = 2.5 -> maybe less clustering ?
%simName = 'VELTEST-2013-12-28-12-53/'; % like above, initV = 5 -> even here they cluster! A lot of bouncing, but then it's ok
% simName = 'VELTEST-2013-12-28-12-55/'; % like above, initV = 10 -> also here!
% simName = 'VELTEST-2013-12-28-12-57/'; % like above with R = 0.2 and initV = 10 -> they have one big cluster, but far away form truth
% simName = 'VELTEST-2013-12-28-13-0/'; % like above with R = 0.2 and initV = 20 -> still one big cluster far away from truth
% simName = 'VELTEST-2013-12-28-13-2/'; % like above with R = 0.2 and initV = 100!! -> 3 big clusters split up, one in the middle and two at the sides
% % simName = 'VELTEST-2013-12-28-13-4/'; % like above with R = 0.2 and initV = 1000!! -> 3 big clusters
% simName = 'VELTEST-2013-12-28-13-10/'; % like above with R = 0.2 and initV = 10000!! -> a lot of movements, but coordinated. In the end somehow only cluster forms.
% % simName = 'VELTEST-2013-12-28-13-17/'; % exactly like above -> one big string of scientists aligned vertically or horizontally. They go up or down so quickly that they stay in the middle.
%simName = 'VELTEST-2013-12-28-13-26/'; % again velocity 100 -> hit the borders and cluster
simName = 'VELTEST-2013-12-28-13-29/'; % again velocity 100, R=0.9 -> one cluster somewhere
simName = 'VELTEST-2013-12-28-13-31/'; % again velocity 100, R=0.02 -> lot of chaos at the beginning, but then slowely clusters forms

% simName = 'SIGMATEST-2014-1-9-14-0/'; % sigma = 0.1, higher than usual. Still cluster forms (R = 0.02, alpha = 0.1)
% simName = 'SIGMATEST-2014-1-9-14-3/'; % sigma = 0.2, higher than usual. Still cluster forms, and they converge to the center very quickly. Velocity at the end is very low.
% simName = 'SIGMATEST-2014-1-9-14-5/'; % sigma = 0.3, higher than usual. Still cluster forms, and they converge to the center very quickly. Velocity at the end is very low.
% simName = 'SIGMATEST-2014-1-9-14-7/'; % sigma = 0.4, higher than usual. Still cluster forms, and they converge to the center very quickly. Velocity at the end is very low.
% simName = 'SIGMATEST-2014-1-9-14-10/'; % sigma = 0.5, higher than usual. Still cluster forms, and they converge to the center very quickly. Velocity at the end is very low.
% simName = 'SIGMATEST-2014-1-9-14-12/'; % sigma = 1, higher than usual. Still cluster forms, and they converge to the center very quickly. Velocity at the end is very low.
% simName = 'SIGMATEST-2014-1-9-14-14/'; % sigma = 1, AND v = 10 higher than usual. Still cluster forms, and they converge to the center very quickly. Velocity at the end is very low.
% simName = 'SIGMATEST-2014-1-9-14-18/'; % sigma = 1, AND epsilon = 0.3. Still cluster forms, and they converge to the center very quickly. Velocity at the end is very low.
% simName = 'SIGMATEST-2014-1-9-14-19/'; % sigma = 1, AND epsilon = 0.8. Still cluster forms, and they converge to the center very quickly. Velocity at the end is very low.
% simName = 'SIGMATEST-2014-1-9-14-22/'; % sigma = 1, AND epsilon = 0.8 AND R = 0. Converge to the center very quickly. Velocity constant.
% simName = 'SIGMATEST-2014-1-9-14-26/'; % sigma = 1, AND epsilon = 0.8 AND R = 0.5. Converge to the center very quickly. Velocity constant.

% Conclusion: increasing sigma is always beneficial for finding the truth.


% simName = 'INDI-2014-1-9-14-29/'; % (R = 0.02, alpha = 0.1, epsilon = 0.1, sigma = 0.01) clustering forming
% simName = 'INDI-2014-1-9-14-31/'; % not_indi, but normal, as above. no main differences, it seems.

simName = 'BIGR-2014-1-23-10-5/'; % R = 0.3, alpha 0.5 v = 1
simName = 'BIGR-2014-1-23-10-7/'; % v = 10
simName = 'BIGR-2014-1-23-10-10/'; % alpha = 0.01
simName = 'BIGR-2014-1-23-10-12/'; % R = 1
simName = 'BIGR-2014-1-23-10-15/'; % v = 0.2
simName = 'BIGR-2014-1-23-10-29/'; % sigma = 0
simName = 'BIGR-2014-1-23-10-30/'; % v = 10
simName = 'BIGR-2014-1-23-10-47/'; % v = 10


% sigma = 0.1, alpha = 0.01, vscaling = 10;
simName = 'SPEEDTEST-2014-1-23-11-23/'; % R = 1 -> 355
simName = 'SPEEDTEST-2014-1-23-11-24/'; % R = 0.02 -> 583

%DUMPDIR = '/mnt/tmp/dump/NAVNP/';


% CLUSTERs vs CLUSTER
simName = 'CLA-2014-2-8-13-43/'; % just checking if 4 clusters are faster than 1 cluster
simName = 'CLA-2014-2-8-13-46/'; % just checking if 4 clusters are faster than 1 cluster
simName = 'CLA-2014-2-8-13-58/'; % 10 random clusters
simName = 'CLA-2014-2-8-14-31-a/'; % 10 and 11 random clusters
simName = 'CLA-2014-2-8-15-7/'; % 20 and 10 clusters on a circle
simName = 'CLA-2014-2-8-15-12/'; % 20 and 10 clusters on a circle, R = 0.03
simName = 'CLA-2014-2-8-15-24/'; % 20 and 10 clusters on a circle, R = 0.03, TAU = 1
%simName = 'CLA-2014-2-8-15-35/'; % 1 cluster on a circle, R = 0.03
simName = 'CLA-2014-2-8-15-59-a/'; % 20 and 10 clusters on a circle, R = 0.03, TAU = 1
simName = 'CLA-2014-2-8-16-6/'; % 20 and 10 clusters on a circle, R = 0.03, TAU = 1

simName = 'LOWALPHA-2014-2-8-18-41/'; % R = 0.03, alpha = 0.01
simName = 'LOWALPHA-2014-2-8-18-46/'; % R = 0.03, alpha = 0.99


% LINEAR
%simName = 'attrLinear_navnp_RClean_n100_fv0_s1_epsilon/attrLinear_navnp_RClean_n100_fv0_s1_epsilon_s0/';
%simName = 'attrLinear_navnp_RClean_n100_fv0_s1_epsilon_v/attrLinear_navnp_RClean_n100_fv0_s1_epsilon_v_s1/';


% MILLEAN
% simName = 'attrMillean_navnp_RClean_n100_fv0_s1_epsilon_v/attrMillean_navnp_RClean_n100_fv0_s1_epsilon_v_s0/';


% ZERO
% simName = 'attrZero_navnp_RClean_n100_fv0_s1_epsilon_v/attrZero_navnp_RClean_n100_fv0_s1_epsilon_v_s1/';

% HARD
%simName = 'attrHard_navnp_RClean_n100_fv0_s1_epsilon/attrHard_navnp_RClean_n100_fv0_s1_epsilon_s5/';


simName = 'AAA-2014-2-7-9-39/';


simName = 'SCENARIO-2014-2-9-16-52/'; % few clusters converge in one
simName = 'SCENARIO-2014-2-9-17-21/'; % alpha = .98 R = 0.03, sigma = 0.02, epsilon = 0.25. Few clusters in a noisy landscape
simName = 'SCENARIO-2014-2-9-17-27/'; % alpha = .99 R = 0.03, sigma = 0.01, epsilon = 0.1. Many singles and many small clusters

simName = '/LOW_TRUTH-2014-2-10-14-53/';
simName = '/LOW_TRUTH-2014-2-10-14-54/';

simName = '/NEWSIMA/';

% TAU test
simName = 'AAA-2014-2-7-9-20/'; % t = 1 cluster near truth
simName = 'AAA-2014-2-7-9-23/'; % t = 2 cluster near truth (a bit further)
simName = 'AAA-2014-2-7-9-24/'; % t = 3 slower convergence, going a bit further
simName = 'AAA-2014-2-7-9-26/'; % t = 4 slower convergence, do not see how far it goes
simName = 'AAA-2014-2-7-9-28/'; % (t_end = 20) t = 4 slower convergence, it stops quite far away
simName = 'AAA-2014-2-7-9-34/'; % (t_end = 20) t = 10 really far away
simName = 'AAA-2014-2-7-9-39/'; % (t_end = 40) t = 10 really far away

simName = 'TTAU-2014-2-20-9-3/'; % t = 13
simName = 'TTAU-2014-2-20-14-15/'; 

% DUMPDIR = '/home/stefano/HS/';
% simName = 'final_tau_20000/final_tau_20000_s1/';


simName = 'EEH-2014-1-27-20-41/'; % v = 1, tau = 10, agents go very slowly to Truth, they don't have enough time
simName = 'EEH-2014-1-27-20-43/'; % v = 1, tau = 1000, agents do not even move...
% simName = 'EEH-2014-1-27-20-46/'; % v = 1, tau = 0.01, agents immediately on truth
% simName = 'EEH-2014-1-27-20-49/'; % v = 1, tau = 0.1, agents go quickly on truth
% simName = 'EEH-2014-1-27-20-57/'; % v = 100, tau = 0.01.
% simName = 'EEH-2014-1-27-20-59/'; % v = 1000, tau = 0.01.
% simName = 'EEH-2014-1-27-21-5/'; % v =1 tau = 2, R = 0.03. clusters forming far away, going to truth
% simName = 'EEH-2014-1-27-21-8/'; % v = 1 tau = 2, R = 0.6. slowly to truth
% simName = 'EEH-2014-1-27-21-9/'; % v = 3 tau = 2, R = 0.6. 1 cluster forming a bit far away, going into truth
% simName = 'EEH-2014-1-27-21-11/'; % v = 5 tau = 2, R = 0.6. too slow..
simName = 'EEH-2014-1-27-21-18/'; % v = 5 tau = 2, R = 0.6. alpha = .99
simName = 'EEH-2014-1-27-21-23-a/'; % v = 10 tau = 1, R = 0.6. alpha = .99 chaos, when agents slow down, 1 cluster on truth
simName = 'EEH-2014-1-27-21-25/'; % v = 10 tau = 1, R = 0.6. alpha = 0.2, 1 cluster mid-far away, then it goes closer
simName = 'EEH-2014-1-27-21-26/'; % v = 100 tau = 1, R = 0.6. alpha = 0.2
simName = 'EEH-2014-1-27-21-28/'; % v = 100 tau = 1, R = 0.6. alpha = 0.2

simName = 'EEH-2014-1-27-22-13/'; % v = 100 tau = 1, R = 0.03. alpha = 0.2, sigma = 0
simName = 'EEH-2014-1-27-22-15/'; % v = 1 tau = 1, R = 0.3. alpha = 0.2, sigma = 0
simName = 'EEH-2014-1-27-22-32/'; % 
simName = 'EEH-2014-1-27-22-34/'; % 
simName = 'EEH-2014-1-27-22-36/'; % 
simName = 'EEH-2014-1-27-22-38/'; % 
simName = 'EEH-2014-1-27-22-40/'; % 
simName = 'EEH-2014-1-27-22-45/'; % 

simName = 'EEH-2014-1-28-0-11/'; % 
simName = 'EEH-2014-1-28-0-15/'; % 





%RBAND
simName = 'RBAND-2014-3-24-21-58/'; % R = 0.3
simName = 'RBAND-2014-3-25-10-36/'; % R = 0.03

simName = 'RBAND-2014-3-25-10-45/'; % R = 0.03, alpha = 0.99
simName = 'RBAND-2014-3-25-10-49/';


%TRUTH-ASIDE

%DUMPDIR = '/home/stefano/Documents/mypapers/swarm_science/data/truth_aside_R_alpha/'
simName = 'truth_aside_R_alpha_s1/';

% Other

simName = 'AAA-2014-2-7-9-39/';

simName = 'NAVPEPSILON-2013-12-21-13-19/';

simName = 'SNAP-2014-5-4-23-12/'; % many clusters far away from truth, boundary effects
simName = 'SNAP-2014-5-4-23-22/'; % many clusters, not so beautiful though
simName = 'SNAP-2014-5-4-23-28/'; % as before

simName = 'SNAP-2014-5-5-9-38/'; % 

simName = 'BIGR-2014-1-23-10-5/';

% USHAPE

simName = 'USHAPE-2014-5-7-15-29/'; % alpha = 0.01
simName = 'USHAPE-2014-5-7-15-39/'; % alpha = 0.5


simName = 'USHAPE-2014-5-7-15-33/'; % alpha = 0.000001

simName = 'EEH-2014-1-27-20-43/'; % v = 1, tau = 1000, agents do not even move...

simName = 'EEH-2014-1-27-20-41/'; % v = 1, tau = 10, agents go very slowly to Truth, they don't have enough time

% TAU test

simName = 'AAA-2014-2-7-9-26/'; % t = 4 slower convergence, do not see how far it goes

simName = 'AAA-2014-2-7-9-39/'; % (t_end = 40) t = 10 really far away



simName = 'TTAU-2014-3-12-11-15/' % t=10, but small R (0.03)

simName = 'TTAU-2014-3-12-11-30/' % t=50, but small R (0.03)


simName = 'TTAU-2014-3-12-14-19/' % t=10, big R (0.3), time2complete 150!


simName = 'TTAU-2014-3-12-13-3/' % t=50, big R (0.3), time2complete 150!


simName = 'AAA-2014-2-7-9-20/'; % t = 1 cluster near truth

simName = 'AAA-2014-2-7-9-23/'; % t = 2 cluster near truth (a bit further)


simName = 'AAA-2014-2-7-9-24/'; % t = 3 slower convergence, going a bit further

simName = 'AAA-2014-2-7-9-28/'; % (t_end = 20) t = 4 slower convergence, it stops quite far away
simName = 'AAA-2014-2-7-9-34/'; % (t_end = 20) t = 10 really far away



simName = 'TTAU-2014-2-20-14-15/'; % t=10, but small R (0.03)

simName = 'TTAU-2014-2-20-9-3/'; % t = 13

simName = 'TTAU-2014-3-12-12-5/' % t=50, big R (0.3)




simName = 'TAUEFFECT-2014-5-15-16-4/' % t = 50, small 5, big alpha


simName = 'TAUEFFECT-2014-5-15-16-0/' % t = 50, small 5, small alpha


simName = 'TAUEFFECT-2014-5-15-17-6/'; % with no noise, no boundaries

simName = 'TAUEFFECT-2014-5-16-9-20/'; % t = 100, big R small alpha, no boundaries

simName = 'TAUEFFECT-2014-5-15-16-44/'; % t = 100, big R small alpha

dumpDir = [DUMPDIR simName];
%myFile = '1291-1';

myFile = '1-1'; % 1521

if (SINGLE)
    videoFile = [VIDEODIR videoSubDir myFile '.avi'];
    if (SAVE_VIDEO && exist([VIDEODIR videoSubDir],'dir') == 0)
        mkdir([VIDEODIR videoSubDir]);
    end
    makevideo(dumpDir, myFile, SAVE_VIDEO, videoFile, plottype, SHOW_POTENTIAL);
    return
end

myFiles = { ...
   '1-1.mat', '3-1.mat', '7-1.mat', '30-1.mat', '2970-1.mat' ...
};

for i=1:length(myFiles)  
    videoFile = [VIDEODIR videoSubDir myFiles{i} '.avi'];
    myFiles{i}
    if (SAVE_VIDEO && exist([VIDEODIR videoSubDir],'dir') == 0)
        mkdir([VIDEODIR videoSubDir]);
    end

    makevideo(dumpDir, myFiles{i}, SAVE_VIDEO, videoFile, plottype, SHOW_POTENTIAL);
end

toc