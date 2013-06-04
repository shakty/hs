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
PLOTS = 1;

DUMPDIR = 'dump/';

simName = 'test_t-2013-6-4-12-1/';
 
temporal_analysis(DUMPDIR, simName, PRECISION, CSV_DUMP, PLOTS);

%%




