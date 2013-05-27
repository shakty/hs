%% Perform Cluster Analysis on the dump

close all
clc
clear

path(path,'util/');

VIDEO = 1;
MPEG  = 0;
DEBUG = 0;

DUMPDIR = '/mnt/tmp/dump/';


%'the_loop-2013-3-8-16-13/';
%MYDIR = 'circle_maybe-2013-3-8-13-22/';
%MYDIR = 'many_little_clusters_do_not_find_truth_mod_a_little_bigger_clusters-2013-3-8-17-29/'
%MYDIR = 'few_big_groups_do_not_find_truth_in_a_smaller_space_cluster_together-2013-3-8-18-28/';

%MYDIR = 'few_big_groups-DIM-vs-ALPHA/';

% WITH THE BUG (full) 'the_loop-2013-3-8-16-13/'
% WITH THE BUG (withtout A and B) 'the_loop-2013-3-9-16-23/'
% WITHOUT THE BUG (no alpha - R): 'the_loop-2013-3-9-16-29/'
% the_loop-2013-3-9-16-42-a/


MYDIR = 'R-alpha-noA-noB-2013-3-6-20-8/'; 

MYDIR = 'tests/the_loop-2013-3-9-22-40/';

MYDIR = 'tests/motus_perpetuous-2013-3-11-12-40-a/';

MYDIR = 'tests/quick_convergence-2013-3-11-12-54/';

MYDIR = 'tests/average_clusters-2013-3-11-13-2/';

MYDIR = 'tests/big_clusters-2013-3-11-13-8/';

MYDIR = 'tests/slow_formation_of_clusters-2013-3-11-14-3/';

MYDIR = 'limited_sigma_R_seq_rnd_avv1/';

MYDIR = 'alpha1_tau_vinit_av0/';

MYDIR = 'cluster_zone_sigma_R_alpha/';

MYDIR = 'truth_corner_alpha_R_av1/';


MYDIR = 'simul_tcorner_noisev_av1_rleft/';

dumpDir = [DUMPDIR MYDIR]

VIDEODIR = '/home/stefano/hs/videos/';


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




videoSubDir = 'simul_tcorner_noisev_av1_rleft/';
videoFile = [VIDEODIR videoSubDir 'few_little_cluster_R=.07_n=0_a=.4.avi'];
makevideo([dumpDir '1592-1.mat'], 0, videoFile);


ded = ded

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
