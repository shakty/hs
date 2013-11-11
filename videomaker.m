%% Perform Cluster Analysis on the dump

close all
clc
clear

path(path,'util/');


DUMPDIR = '/mnt/tmp/dump/NEW/';
MYDIR = 'attrZero_nv_Kseed_rndseq_tm_Rleft';
MYDIR = 'attrLinear_nv_rndseq_tm_Rleft';

SIGMA = '_s5';

dumpDir = [DUMPDIR MYDIR '/' MYDIR SIGMA '/'];

VIDEODIR = '/home/stefano/hs/videos/NEW/';


SAVE_VIDEO = 0;
SHOW_POTENTIAL = 1;
plottype = 0;

% PLOT TYPE
% plot_cross = 0;
% plot_number = 1;
% plot_number_color = 2;
% plot_arrow = 3;
    
myFile = '3131-1.mat';
videoSubDir = [ MYDIR '/' ];
videoFile = [VIDEODIR videoSubDir myFile '.avi'];
if (exist([VIDEODIR videoSubDir],'dir') == 0)
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
