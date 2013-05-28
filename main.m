%% Stefano Balietti

%% Initialization of variables 
%Clear workspace
close all
clear
clc

%% Set Random Seed

%s = RandStream('mcg16807','Seed',0);
%RandStream.setGlobalStream(s)

%% Add other directories to path
path(path,'util/'); % Help functions
path(path,'lib/'); % Help functions

% Default Values

simName = 'Sim';
dumpDir = 'dump/';
confDir = 'conf/';

compLOCAL = 0;
compPARALLEL = 1;
compLSF = 2;
COMPUTATION = compLOCAL;

VIDEO = 0;
DEBUG = 1;
DUMP = 1;

% CHANGE AFTER!!!!!

VIDEO = 1;
DUMP = 0;
COMPUTATION = compLOCAL;

% Force Local Computation
%COMPUTATION = compLOCAL;

% Force Local Computation

nRuns = 1;
t_ends = 1;
As = [0];
Bs = [0];

%n_agents  = 2;

taus = 2;
vScalings = [1];
alphas = [0.05];       	% weighting of velocity terms
Rs     = [0.03];
truths = [0.5;0.5];
ideas_space_sizes = 1;
ks=1;
sigmas = 0.1;
n_agents = 100;

load([confDir 'refactor']);
%simParamsStruct.VIDEO=1;
simParamsStruct.simName = createSimName(simParamsStruct.simName,simParamsStruct.DUMP,simParamsStruct.dumpDir);

%set(gcf, 'DoubleBuffer', 'on');

%% Start Vectorization of Parameters Sets

switch (COMPUTATION)
    case compLSF
         param_sets_LSF (simParamsStruct);
        
%    case compPARALLEL
%        param_sets_parallel (dumpDir,simName,VIDEO,DEBUG,DUMP,...
%            nRuns, dts,t_ends,n_agents,ideas_space_sizes,ideas_space_dims,...
%            As,Bs,ks,d0s,d1s,alphas,taus,Rs,sigmas,...
%            vScalings,nClusters,clusterTightness,truths);
        
    case compLOCAL
        param_sets_local (simParamsStruct);
end

fprintf('\n%s: execution completed correctly\n', simParamsStruct.simName);
% Exit Matlab when invoked from command line with -r option
%exit
