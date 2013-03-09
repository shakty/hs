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
DIR = 'circle_maybe-2013-3-8-13-22';
DIR = 'tests/the_loop-2013-3-9-22-40/'

tocsv(DIR, CSV_CLU, CSV_POS, PLOT_POS, PLOT_CLU);

%%




