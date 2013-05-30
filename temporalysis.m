%% Saves simulations into properly formatted CSV files

close all;
clear;
clc;

%% Add other directories to path
path(path,'util/'); % Help functions
path(path,'lib/'); % Help functions


PRECISION = 100;

CSV_DUMP = 1;
PLOTS = 0;

DUMPDIR = 'dump/refactor/';
simName = 'refactor-2013-5-29-14-52';
 
temporal_analysis(DUMPDIR, simName, PRECISION, CSV_DUMP, PLOTS);

%%




