%% Stefano Balietti

%% Initialization of variables 
%Clear workspace
close all
clear
clc

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

VIDEO = 1;
DEBUG = 0;
DUMP = 1;

% Conf to execute
%movingTau_LOCAL;
%path(path,'conf/'); % Configuration files 
%velocity_still_LOCAL;

%load([confDir 'few_big_groups_do_not_find_truth']);

MYSIM = 'size-alpha-R-tau-vscaling';
MYSIM = 'sigma_tau';

load([confDir 'OK/' MYSIM]);

%load([confDir 'aa']);


% Force Local Computation
 COMPUTATION = compLOCAL;

simName = createSimName(simName,DUMP,dumpDir);

%set(gcf, 'DoubleBuffer', 'on');

%% Start Vectorization of Parameters Sets

switch (COMPUTATION)
    case compLSF
         param_sets_LSF (dumpDir,simName,VIDEO,DEBUG,DUMP,...
            nRuns, dts,t_ends,n_agents,ideas_space_sizes,ideas_space_dims,...
            As,Bs,ks,d0s,d1s,alphas,taus,Rs,sigmas,...
            vScalings,nClusters,clusterTightness,truths);
        
    case compPARALLEL
        param_sets_parallel (dumpDir,simName,VIDEO,DEBUG,DUMP,...
            nRuns, dts,t_ends,n_agents,ideas_space_sizes,ideas_space_dims,...
            As,Bs,ks,d0s,d1s,alphas,taus,Rs,sigmas,...
            vScalings,nClusters,clusterTightness,truths);
        
    case compLOCAL
        param_sets_local (dumpDir,simName,VIDEO,DEBUG,DUMP,...
            nRuns, dts,t_ends,n_agents,ideas_space_sizes,ideas_space_dims,...
            As,Bs,ks,d0s,d1s,alphas,taus,Rs,sigmas,...
            vScalings,nClusters,clusterTightness,truths);
end

fprintf('\n%s: execution completed correctly\n',simName);
% Exit Matlab when invoked from command line with -r option
%exit
