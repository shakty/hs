%% Perform Cluster Analysis on the dump

close all
clc
clear

path(path,'util/');

VIDEO = 1;
MPEG  = 1;
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

MYDIR = 'cluster_zone_sigma_R_alpha/';

MYDIR = 'alpha1_tau_vinit_av1/';

dumpDir = [DUMPDIR MYDIR]

VIDEODIR = '/tmp/';%dumpDir;

simNumber = 495;
simCount = 1;

makevideo([dumpDir '13-1.mat'],0,0)


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
    makevideo(fileName, fileOut, MPEG);
end
