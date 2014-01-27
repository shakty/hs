%% Aggregates the results of the analysis of the simulation results.
tic;

close all;
clear;
clc;

%% Add other directories to path
path(path,'util/'); % Help functions
path(path,'lib/'); % Help functions

DUMPDIR = '/home/stefano/HS/attrLinear_navnp_SpeedTest_fullscan/';

DUMPDIR = '/home/stefano/hs/test/EEH-2014-1-27-22-58/';

speedtest_fun(DUMPDIR);