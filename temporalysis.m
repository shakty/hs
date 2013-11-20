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


% Dump only every X steps
DUMP_RATE = 100;

% Only clusters of size above the cutoff are included in the analysis
CLU_CUTOFF = 2;
% When computing the coverage we build a grid on top of the space of cell
% size = PRECISION
PRECISION = 100;

CSV_DUMP = 0;
PLOTS = 0;

DUMPDIR = '/mnt/tmp/dump/RND_SEED/';
simName = 'attrZero_nav_rndseeds_rndseq_tm_Alpha1_n100_fv0/attrZero_nav_rndseeds_rndseq_tm_Alpha1_n100_fv0_s1';

DUMPDIR = '/tmp/';
simName = 'aa';

DUMPDIR = '/home/stefano/hs/test/';
simName = 'TEST_SIZE10-2013-11-20-22-10/';

dumpDir = [DUMPDIR simName '/']

tic
temporal_analysis(DUMPDIR, simName, PRECISION, CLU_CUTOFF, CSV_DUMP, DUMP_RATE, PLOTS);
toc
