%% Saves simulations into properly formatted CSV files

close all;
clear;
clc;

%% Add other directories to path
path(path,'util/'); % Help functions
path(path,'lib/'); % Help functions

% Change default axes fonts.
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 14)

% Change default text fonts.
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 14)

PRECISION = 100;

CSV_DUMP = 1;
PLOTS = 0;

DUMPDIR = 'dump/refactor/';
simName = 'refactor-2013-5-31-12-26';
 
temporal_analysis(DUMPDIR, simName, PRECISION, CSV_DUMP, PLOTS);

%%




