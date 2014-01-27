%% Aggregates the results of the analysis of the simulation results.
tic;

close all;
clear;
clc;

%% Add other directories to path
path(path,'util/'); % Help functions
path(path,'lib/'); % Help functions

DUMPDIR = '/home/stefano/HS/attrLinear_navnp_SpeedTest_fullscan/';

speedtest_fun(DUMPDIR);