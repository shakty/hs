%% Saves simulations into properly formatted CSV files

close all;
clear;
clc;

%% Add other directories to path
path(path,'../'); % Help functions
path(path,'../util/'); % Help functions
path(path,'../lib/'); % Help functions

% Change default axes fonts.
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 14)

% Change default text fonts.
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 14)

% Dump to CSV only every X steps. 
DUMP_RATE = 1; % 100 should produces ~22 iterations

PRECISION = 100;

DUMP = 0;
PLOTS = 0;

DUMPDIR = '/home/stefano/hs/test/NEWTEST-2013-12-8-17-49/';
simName = 'A/';

dumpDir = [DUMPDIR simName];

outDir = [dumpDir '/' 'agents/'];
% Creating outDir if not existing.
if (exist(outDir, 'dir') == 0)
    mkdir(outDir);
end

fileName = '1-1.mat';

paramsObj = struct( ...
                    'folderName', DUMPDIR, ...
                    'simName', simName, ...
                    'fileName', fileName, ...
                    'outDirAgents', outDir, ...
                    'DUMP', DUMP, ...
                    'DUMP_RATE', DUMP_RATE, ...
                    'PLOTS', PLOTS, ...
                    'PRECISION', PRECISION ...
);

tic
agents_onefile(paramsObj);
toc