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
plottype = 0;
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
simName = 'VELTEST-2013-12-28-12-43/'; % like above, initV = 0.2 -> 1 big cluster and a few outsiders
% simName = 'VELTEST-2013-12-28-12-45/'; % like above, initV = 0.1 -> 1 big cluster
% simName = 'VELTEST-2013-12-28-12-47/'; % like above, initV = 1.5 -> many small clusters, some bouncing
% simName = 'VELTEST-2013-12-28-12-49/'; % like above, initV = 2 -> more bouncing, many more spread out clusters
% simName = 'VELTEST-2013-12-28-12-50/'; % like above, initV = 2.5 -> maybe less clustering ?
% simName = 'VELTEST-2013-12-28-12-53/'; % like above, initV = 5 -> even here they cluster! A lot of bouncing, but then it's ok
% simName = 'VELTEST-2013-12-28-12-55/'; % like above, initV = 10 -> also here!
% simName = 'VELTEST-2013-12-28-12-57/'; % like above with R = 0.2 and initV = 10 -> they have one big cluster, but far away form truth
% simName = 'VELTEST-2013-12-28-13-0/'; % like above with R = 0.2 and initV = 20 -> still one big cluster far away from truth
% simName = 'VELTEST-2013-12-28-13-2/'; % like above with R = 0.2 and initV = 100!! -> 3 big clusters split up, one in the middle and two at the sides
% % simName = 'VELTEST-2013-12-28-13-4/'; % like above with R = 0.2 and initV = 1000!! -> 3 big clusters
% % simName = 'VELTEST-2013-12-28-13-10/'; % like above with R = 0.2 and initV = 10000!! -> a lot of movements, but coordinated. In the end somehow only cluster forms.
% % simName = 'VELTEST-2013-12-28-13-17/'; % exactly like above -> one big string of scientists aligned vertically or horizontally. They go up or down so quickly that they stay in the middle.
% simName = 'VELTEST-2013-12-28-13-26/'; % again velocity 100 -> hit the borders and cluster
% simName = 'VELTEST-2013-12-28-13-29/'; % again velocity 100, R=0.9 -> one cluster somewhere
% simName = 'VELTEST-2013-12-28-13-31/'; % again velocity 100, R=0.02 -> lot of chaos at the beginning, but then slowely clusters forms


DUMPDIR = '/mnt/tmp/dump/NAVNP/';
simName = 'attrLinear_navnp_RClean_n100_fv0_s1_epsilon/attrLinear_navnp_RClean_n100_fv0_s1_epsilon_s1/';

dumpDir = [DUMPDIR simName];
%myFile = '1291-1';

myFile = '21-1';

if (SINGLE)
    videoFile = [VIDEODIR videoSubDir myFile '.avi'];
    if (SAVE_VIDEO && exist([VIDEODIR videoSubDir],'dir') == 0)
        mkdir([VIDEODIR videoSubDir]);
    end
    makevideo(dumpDir, myFile, SAVE_VIDEO, videoFile, plottype, SHOW_POTENTIAL);
    return
end

myFiles = {
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