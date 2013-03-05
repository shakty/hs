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

DIR = 'testcsv-2013-3-5-9-47';

tocsv(DIR, CSV_CLU, CSV_POS, PLOT_POS, PLOT_CLU);

%%




