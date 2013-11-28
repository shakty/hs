%% Perform Cluster Analysis on the dump

close all
clc
clear

path(path,'util/');

% EXP_BUG

DUMPDIR = '/mnt/tmp/dump/EXP_BUG/';

% R CLEAN
MYDIR = 'attrHard_nav_rndseeds_rndseq_tm_RClean_n100_fv0';
MYDIR = 'attrZero_nav_rndseeds_rndseq_tm_RleftClean_n100_fv0';
MYDIR = 'attrLinear_nav_rndseeds_rndseq_tm_Rclean_n100_fv0';
MYDIR = 'attrMillean_nav_rndseeds_rndseq_tm_RClean_n100_fv0';

% R = 0;
MYDIR = 'attrLinear_nav_rndseeds_rndseq_tm_R0_n100_fv0';
%MYDIR = 'attrLinear_nav_rndseeds_rndseq_tm_Alpha1_n100_fv0';
MYDIR = 'attrZero_nav_rndseeds_rndseq_tm_R0_n100_fv0';
MYDIR = 'attrHard_nav_rndseeds_rndseq_tm_R0_n100_fv0';
MYDIR = 'attrMillean_nav_rndseeds_rndseq_tm_R0_n100_fv0';

% SIZE 10

% R CLEAN
MYDIR = 'attrZero_nav_rndseeds_rndseq_tm_RClean_n100_fv0_s10';
MYDIR = 'attrLinear_nav_rndseeds_rndseq_tm_RClean_n100_fv0_s10';
MYDIR = 'attrMillean_nav_rndseeds_rndseq_tm_RClean_n100_fv0_s10';
% MYDIR = 'attrHard_nav_rndseeds_rndseq_tm_RClean_n100_fv0_s10'; (TODO);

DUMPDIR = '/mnt/tmp/dump/SIZE10/';


SIGMA = '_s0';

% DUMP DIR
dumpDir = [DUMPDIR MYDIR '/' MYDIR SIGMA '/'];
% FILE
myFile = '1000-1.mat';


VIDEODIR = '/home/stefano/hs/videos/NEW/';
SAVE_VIDEO = 0;
SHOW_POTENTIAL = 0;
plottype = 0 ;

% PLOT TYPE
% plot_cross = 0;
% plot_number = 1;
% plot_number_color = 2;
% plot_arrow = 3;
    

videoSubDir = [ MYDIR '/' ];
videoFile = [VIDEODIR videoSubDir myFile '.avi'];
if (SAVE_VIDEO && exist([VIDEODIR videoSubDir],'dir') == 0)
    mkdir([VIDEODIR videoSubDir]);
end
makevideo([dumpDir myFile], SAVE_VIDEO, videoFile, plottype, SHOW_POTENTIAL);

return;

% parentDir = '/home/stefano/hs/dump/tests/';
%myFiles = {
%    'cluster_zone_R0.2_av1-2013-3-28-11-33', ...
%    'cluster_zone_R0.2_sigma.1_av0-2013-3-28-12-6', ...
%    'cluster_zone_R0.2_sigma.5_av0-2013-3-28-11-53', ...
%    'cluster_zone_R0.2_sigma.5_av1-2013-3-28-11-40', ...
%    'cluster_zone_R0.2_sigma.9_av0-2013-3-28-11-50', ...
%    'cluster_zone_R0.2_sigma.9_av1-2013-3-28-11-44'
%};
%videoSubDir = 'truth_corner_av0/';
%for i=1:length(myFiles)  
%    videoFile = [VIDEODIR videoSubDir myFiles{i} '.avi'];
%    makevideo([parentDir myFiles{i} '/1-1.mat'], 1, videoFile);
%end

files = dir(dumpDir);
fileIndex = find(~[files.isdir]);



%%

display(['Files Found: ' num2str(length(fileIndex))]);

if (~length(fileIndex) > 1)
    display('Empty dir.');
end

%%

for i = 1:length(fileIndex)
     
    append = files(fileIndex(i)).name;
    fileName = [dumpDir, append];
    
    % We load only .mat
    [PATH,NAME,EXT] = fileparts(fileName);
    if (~strcmpi(EXT,'.mat')) 
        continue;
    end
    
    fileOut = [VIDEODIR append '.avi'];
    makevideo(fileName, MPEG, fileOut);
end
