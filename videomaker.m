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
plottype = 0 ;
SINGLE = 1;

% PLOT TYPE
% plot_cross = 0;
% plot_number = 1;
% plot_number_color = 2;
% plot_arrow = 3;

% SINGLE
DUMPDIR = '/home/stefano/hs/test/'; 
simName = 'TESTLINEAR/';
dumpDir = [DUMPDIR simName];
myFile = '1291-1';


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