%% Stefano Balietti

%% Initialization of variables 
%Clear workspace
close all
clear
%clc

%% Add other directories to path
path(path,'util/'); % Help functions
path(path,'lib/'); % Help functions



    headers_clusters = {
        'sim', ...
        'run', ...
        't', ...
        'count', ...
        'size.avg', ...
        'size.sd', ...
        'fromtruth.avg', ...
        'fromtruth.sd'};

    headers_params = {
        'sim', ...
        'run', ...
        'timestamp', ...
        't.end', ...
        'dt', ...
        'nagents', ...
        'spacesize',...
        'spacedim', ...
        'alpha', ... % Own velocity
        'R', ...
        'k', ... % Exponent forces
        'A', ... % Attractive force
        'd0', ... 
        'B', ... % Repulsive force
        'd1', ... 
        'tau', ... % Truth strength
        'sigma', ... % Std noise
        'init.vscaling', ...
        'init.nclusters', ...
        'init.clusterradio', ...
        'truth.x', ...
        'truth.y', ...
        'conv'};

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

%MYSIM = 'size-alpha-R-tau-vscaling';
%MYSIM = 'sigma_tau';

%load([confDir 'OK/' MYSIM]);

load([confDir 'TESTS/' 'the_loop']);

% FORCE FOR NOW...
VIDEO = 1;
DEBUG = 0;
DUMP = 0;


% Force Local Computation
COMPUTATION = compLOCAL;
nRuns = 1;
t_ends = 10;
As = [1];
Bs = [1]

taus = 2;
vScalings = [0:0.1:3]
alphas = [0.1];       	% weighting of velocity terms
Rs     = [0.1];
truths = [0.5;0.5];
ks=1

dumpDir = 'dump/tests/'
simName = 'xxx';

simName = createSimName(simName,DUMP,dumpDir);

folderName = [dumpDir simName];
[status,message,messageid] = mkdir(folderName);  
%% ToCSV params

   % params (for both)
   paramFileName = [folderName '/params.csv']
   write_csv_headers(paramFileName, headers_params);


% Open CSV cluster

  clusterFileName = [folderName '/clusters.csv']
  write_csv_headers(clusterFileName, headers_clusters);



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
