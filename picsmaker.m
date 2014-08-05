%% Perform Cluster Analysis on the dump
tic
close all
clc
clear

path(path,'util/');


% Defaults for this blog post
width = 3;     % Width in inches
height = 3;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 11;      % Fontsize
lw = 2;      % LineWidth
msz = 8;       % MarkerSize


% Change default axes fonts.
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 30)

% Change default text fonts.
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 30)

% set(gca,'XTick',-3:3); %<- Still need to manually specific tick marks
% set(gca,'YTick',0:10); %<- Still need to manually specific tick marks
% set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

% The properties we've been using in the figures
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz


% SIGMA
SIGMA = '_s1';

% PARAMS

% plot_cross = 0;
% plot_number = 1;
% plot_number_color = 2;
% plot_arrow = 3;


outDir = '/home/stefano/hs_writeup/imgs/new/snaps_test/';

myDumpDirs = {};
myDirs = {};
myFiles = {};

DUMP_DIR = '/mnt/tmp/dump/';

% myDumpDirs{1} = 'EXP_BUG/';
% myFiles{1} = {'1-1.mat', '3-1.mat', '7-1.mat', '30-1.mat', '2970-1.mat'};
% myDirs{1} = {
%     'attrLinear_nav_rndseeds_rndseq_tm_Rclean_n100_fv0', ...
%     'attrMillean_nav_rndseeds_rndseq_tm_RClean_n100_fv0', ...
%     'attrHard_nav_rndseeds_rndseq_tm_RClean_n100_fv0', ...
%     'attrZero_nav_rndseeds_rndseq_tm_RleftClean_n100_fv0' ...
% };
% 
% myDumpDirs{2} = 'EXP_BUG/';
% myFiles{2} = {'1-1.mat', '3-1.mat', '7-1.mat', '30-1.mat', '99-1.mat'};
% myDirs{2} = {   
%     'attrLinear_nav_rndseeds_rndseq_tm_R0_n100_fv0', ...
%     'attrZero_nav_rndseeds_rndseq_tm_R0_n100_fv0', ...
%     'attrHard_nav_rndseeds_rndseq_tm_R0_n100_fv0', ...
%     'attrMillean_nav_rndseeds_rndseq_tm_R0_n100_fv0' ...
% };
% 
% myDumpDirs{3} = 'SIZE10/';
% myFiles{3} = {'1-1.mat', '3-1.mat', '7-1.mat', '30-1.mat', '2970-1.mat'};
% myDirs{3} = {   
%     'attrZero_nav_rndseeds_rndseq_tm_RClean_n100_fv0_s10', ...
%     'attrLinear_nav_rndseeds_rndseq_tm_RClean_n100_fv0_s10', ...
%     'attrMillean_nav_rndseeds_rndseq_tm_RClean_n100_fv0_s10', ...
%     ... % 'attrHard_nav_rndseeds_rndseq_tm_RClean_n100_fv0_s10'; (TODO);
% };


% CONV FAR FROM TRUTH
DUMP_DIR = '/home/stefano/hs/';
myDumpDirs{1} = 'test/';
myFiles{1} = {'1-1.mat'};
myDirs{1} = {'NOB_TAU_ALPHA-2014-8-4-10-43'};
% 
% % CLUSTERS 
% DUMP_DIR = '/home/stefano/hs/';
% myDumpDirs{1} = 'test/';
% myFiles{1} = {'1-1.mat'};
% myDirs{1} = {'NAVPEPSILON-2013-12-21-13-19/'};
% 
% % CONV ON TRUTH
% DUMP_DIR = '/home/stefano/hs/';
% myDumpDirs{1} = 'test/';
% myFiles{1} = {'1-1.mat'};
% myDirs{1} = {'SNAP-2014-5-5-10-39/'}; % BIGR-2014-1-23-10-5 (just 7 steps)
% 
% % A FEW BIG CLUSTERS THAT CONVERGE IN ONE
% DUMP_DIR = '/home/stefano/hs/';
% myDumpDirs{1} = 'test/';
% myFiles{1} = {'1-1.mat'};
% myDirs{1} = {'SCENARIO-2014-2-9-16-52/'};
 
% NON CONVERGENCE, NO CLUSTERS
% DUMP_DIR = '/home/stefano/hs/';
% myDumpDirs{1} = 'test/';
% myFiles{1} = {'1-1.mat'};
% myDirs{1} = {'SCENARIO-2014-2-9-17-27/'};

% A few big medium clusters NO!
% DUMP_DIR = '/home/stefano/hs/';
% myDumpDirs{1} = 'test/';
% myFiles{1} = {'1-1.mat'};
% myDirs{1} = {'SNAP-2014-5-5-9-38/'};


% A few big medium clusters
% DUMP_DIR = '/home/stefano/hs/';
% myDumpDirs{1} = 'test/';
% myFiles{1} = {'1-1.mat'};
% myDirs{1} = {'SNAP-2014-5-4-23-22/'};


WITH_SIGMA = 0;

myFrames = [1, 200, 1000, 2000];

myFrames = [20, 40, 50, 80, 100, 200, 1000, 2000];

LIMITS = 0;

SHOW_TITLE = 0;
SHOW_AXIS = 1;
SHOW_POTENTIAL = 0;
plottype = 0;

for j=1:length(myDirs)
    
    for h=1:length(myDirs{j})
    
        % subDir has no '/' at the end!
        if (WITH_SIGMA)
            subDir = [ myDirs{j}{h} SIGMA];
        else
            subDir = '';
        end
        
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
            makepics(dumpDir, subDir, file, outSubDir, myFrames, plottype, ...
                        SHOW_POTENTIAL, SHOW_TITLE, SHOW_AXIS, LIMITS);
        end
    end
end
toc