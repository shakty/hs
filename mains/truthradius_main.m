%% Saves simulations into properly formatted CSV files

close all;
clear;
clc;

% SEE AGENTS TO FIX.

%% Add other directories to path
path(path,'util/'); % Help functions
path(path,'lib/'); % Help functions

% Change default axes fonts.
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 14)

% Change default text fonts.
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 14)

% Dump to CSV only every X steps. 
DUMP_RATE = 1; % 100 should produces ~22 iterations

% Different radius to try out (around truth).
RADIUSs = [0.01, 0.05, 0.1, 0.25, 0.4];

% An agent must stay for at least STAY_FOR iterations.
STAY_FOR = 1;

% A share of at least CONSENSUS_THRESHOLD agents must stay within the
% smallest radius to be say that there is a consensus on truth.
CONSENSUS_THRESHOLD = 2/3;

% Consenus must on truth must hold for at least CONSENSUS_ON_TRUTH_FOR
% iterations to say that there is really a consensus on truth.
CONSENSUS_ON_TRUTH_FOR = 20;

DUMP = 1;
PLOTS = 0;

DUMPDIR = '/mnt/tmp/dump/EXP_BUG/';
simName = 'attrLinear_nav_rndseeds_rndseq_tm_Rclean_n100_fv0/attrLinear_nav_rndseeds_rndseq_tm_Rclean_n100_fv0_s1';

DUMPDIR = '/home/stefano/hs/test/';
simName = 'NEWTEST-2013-12-8-17-49/';

dumpDir = [DUMPDIR simName];

outDir = [dumpDir '/' 'truthradius/'];
% Creating outDir if not existing.
if (exist(outDir, 'dir') == 0)
    mkdir(outDir);
end

fileName = '1-1.mat';

paramsObj = struct( ...                    
                    'folderName', DUMPDIR, ...
                    'simName', simName, ...
                    'fileName', fileName, ...
                    'RADIUSs', RADIUSs, ...
                    'STAY_FOR', STAY_FOR, ...
                    'CONSENSUS_ON_TRUTH_FOR', CONSENSUS_ON_TRUTH_FOR, ...
                    'CONSENSUS_THRESHOLD', CONSENSUS_THRESHOLD, ...
                    'outDirRadius', outDir, ...
                    'DUMP', DUMP, ...
                    'DUMP_RATE', DUMP_RATE, ...
                    'PLOTS', PLOTS ...
);

tic
truthradius_onefile(paramsObj);
toc