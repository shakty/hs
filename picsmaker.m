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

% SIGMA
SIGMA = '_s1';

% PARAMS

% plot_cross = 0;
% plot_number = 1;
% plot_number_color = 2;
% plot_arrow = 3;

SHOW_POTENTIAL = 0;
plottype = 0;

outDir = '/home/stefano/hs_writeup/imgs/snapshots_nv/';

myDumpDirs = {};
myDirs = {};
myFiles = {};

DUMP_DIR = '/mnt/tmp/dump/';

myDumpDirs{1} = 'EXP_BUG/';
myFiles{1} = {'1-1.mat', '3-1.mat', '7-1.mat', '30-1.mat', '2970-1.mat'};
myDirs{1} = {
    'attrLinear_nav_rndseeds_rndseq_tm_Rclean_n100_fv0', ...
    'attrMillean_nav_rndseeds_rndseq_tm_RClean_n100_fv0', ...
    'attrHard_nav_rndseeds_rndseq_tm_RClean_n100_fv0', ...
    'attrZero_nav_rndseeds_rndseq_tm_RleftClean_n100_fv0' ...
};

myDumpDirs{2} = 'EXP_BUG/';
myFiles{2} = {'1-1.mat', '3-1.mat', '7-1.mat', '30-1.mat', '99-1.mat'};
myDirs{2} = {   
    'attrLinear_nav_rndseeds_rndseq_tm_R0_n100_fv0', ...
    'attrZero_nav_rndseeds_rndseq_tm_R0_n100_fv0', ...
    'attrHard_nav_rndseeds_rndseq_tm_R0_n100_fv0', ...
    'attrMillean_nav_rndseeds_rndseq_tm_R0_n100_fv0' ...
};

myDumpDirs{3} = 'SIZE10/';
myFiles{3} = {'1-1.mat', '3-1.mat', '7-1.mat', '30-1.mat', '2970-1.mat'};
myDirs{3} = {   
    'attrZero_nav_rndseeds_rndseq_tm_RClean_n100_fv0_s10', ...
    'attrLinear_nav_rndseeds_rndseq_tm_RClean_n100_fv0_s10', ...
    'attrMillean_nav_rndseeds_rndseq_tm_RClean_n100_fv0_s10', ...
    ... % 'attrHard_nav_rndseeds_rndseq_tm_RClean_n100_fv0_s10'; (TODO);
};


myFrames = [1, 200, 1000, 2000];

for j=1:length(myDirs)
    
    for h=1:length(myDirs{j})
    
        % subDir has no '/' at the end!
        subDir = [ myDirs{j}{h} SIGMA];
        dumpDir = [ DUMP_DIR myDumpDirs{j} myDirs{j}{h} '/'];
        outSubDir = [ outDir myDumpDirs{j} '/'];
    
        for i=1:length(myFiles{j})
            if (exist(outDir, 'dir') == 0)
                mkdir(outDir);
            end
            if (exist(outSubDir , 'dir') == 0)
                mkdir(outSubDir);
            end
            file = myFiles{j}{i};
            makepics(dumpDir, subDir, file, outSubDir, myFrames, plottype, SHOW_POTENTIAL);
        end
    end
end
toc