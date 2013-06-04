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

% The minimum size of a cluster to be included in the analysis
CLU_CUTOFF = 10;
% When computing the coverage we build a grid on top of the space of cell
% size = PRECISION
PRECISION = 100;

CSV_DUMP = 1;
PLOTS = 0;

DUMPDIR = 'dump/';

simName = 'attrExpo_nv_rndseq_tm_Rleft/'; 

temporal_analysis(DUMPDIR, simName, PRECISION, CLU_CUTOFF, CSV_DUMP, PLOTS);

%%




