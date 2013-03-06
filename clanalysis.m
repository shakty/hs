%% Saves simulations into properly formatted CSV files

close all;
clear;
clc;

%% Add other directories to path
path(path,'util/'); % Help functions
path(path,'lib/'); % Help functions


CSV_CLU = 1;
CSV_POS = 0;

PLOT_POS = 0;
PLOT_CLU= 0;

DIR = 'few_big_groups-DIM-vs-ALPHA';
DIR = 'sigma_tau-2013-3-6-10-36';
tocsv(DIR, CSV_CLU, CSV_POS, PLOT_POS, PLOT_CLU);

%%




