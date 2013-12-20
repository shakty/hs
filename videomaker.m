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

dumpDir = [DUMPDIR simName];
%myFile = '1291-1';

myFile = '1-1';

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