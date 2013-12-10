%% Aggregates the results of the analysis of the simulation results.

close all;
clear;
clc;

%% Add other directories to path
path(path,'util/'); % Help functions
path(path,'lib/'); % Help functions

aggrParams = 1;

DUMPDIR = '/home/stefano/hs/test/';
subDir = 'NEWTEST-2013-12-8-17-49/';

outDir = [DUMPDIR subDir 'aggr/'];


aggregate_agents(DUMPDIR, subDir, outDir, aggrParams)
aggregate_clusters(DUMPDIR, subDir, outDir)