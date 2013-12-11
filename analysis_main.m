%% Aggregates the results of the analysis of the simulation results.
tic;

close all;
clear;
clc;

%% Add other directories to path
path(path,'util/'); % Help functions
path(path,'lib/'); % Help functions

DUMP = 1;
DUMP_RATE = 1; % Dump every x steps

% Only clusters of size above the cutoff are included in the analysis
CLU_CUTOFF = 2;

% When computing the coverage we build a grid on top of the space of cell
% size = PRECISION
PRECISION = 100;

% Different radius to try out (around truth).
RADIUSs = [0.01, 0.05, 0.1, 0.25, 0.4];
nRadiusesPlusOne = length(RADIUSs) + 1;

% An agent must stay for at least STAY_FOR iterations.
STAY_FOR = 1;

% A share of at least CONSENSUS_THRESHOLD agents must stay within the
% smallest radius to be say that there is a consensus on truth.
CONSENSUS_THRESHOLD = 2/3;

% Consenus must on truth must hold for at least CONSENSUS_ON_TRUTH_FOR
% iterations to say that there is really a consensus on truth.
CONSENSUS_ON_TRUTH_FOR = 20;

aggrParams = 1;

DUMPDIR = '/home/stefano/hs/test/';

simName = 'NEWTEST-2013-12-8-17-49/';

path2sim = [DUMPDIR simName];

aggregateSims = 1;

% Every subdirectory of path2sim contains simulations results.
dirs = dir(path2sim);
dirIndex = find([dirs.isdir]);

if (isempty(dirIndex))
    error('Invalid Directory Selected');
end

outDir = [path2sim 'aggr/'];

% Creating outDir if not existing.
if (exist(outDir, 'dir') == 0)
    mkdir(outDir);
end

% Each subdir containing results must be aggregated (they are divided by
% run) and then the aggregated results must be aggregated overall.
for d = 1:length(dirIndex)

    subDir = dirs(dirIndex(d)).name;

    if ( ...
        strcmpi(subDir,'.') ...
        || strcmpi(subDir,'..') ...
        || strcmpi(subDir,'aggr') ...
        || strcmpi(subDir,'img') ...        
    )
        continue;
    end

    subDir = [subDir '/'];
    
    dirPath = [path2sim subDir];   
 
    outDirAgents = [dirPath '/' 'agents/'];
    outDirRadius = [dirPath '/' 'truthradius/'];
    outDirClusters = [dirPath '/' 'clusters/'];

    % Creating outDir if not existing.
    if (exist(outDirAgents, 'dir') == 0)
        mkdir(outDirAgents);
    end

    % Creating outDir if not existing.
    if (exist(outDirRadius, 'dir') == 0)
        mkdir(outDirRadius);
    end

    % Creating outDir if not existing.
    if (exist(outDirClusters, 'dir') == 0)
        mkdir(outDirClusters);
    end
    
    files = dir(dirPath);
    fileIndex = find(~[files.isdir]);

    if (isempty(fileIndex))
        error('Invalid Directory Selected');
    end
    
    % Number of files.
    nFiles = length(fileIndex);
    
    % Load all parameters matrices in one.
    for f = 1:nFiles

        append = files(fileIndex(f)).name;
        fileName = [dirPath append];

        % We load only .mat
        [PATH, NAME, EXT] = fileparts(fileName);
        if (~strcmpi(EXT,'.mat') ||  ~isempty(strfind(NAME, 'sums_')))  
            continue;
        end
    
        params = struct( ...
                'folderName', path2sim, ...
                'simName', subDir, ...
                'fileName', NAME, ...
                'RADIUSs', RADIUSs, ...
                'STAY_FOR', STAY_FOR, ...
                'CONSENSUS_ON_TRUTH_FOR', CONSENSUS_ON_TRUTH_FOR, ...
                'CONSENSUS_THRESHOLD', CONSENSUS_THRESHOLD, ...
                'DUMP', DUMP, ...
                'DUMP_RATE', DUMP_RATE, ...
                'PLOTS', 0, ...
                'CLU_CUTOFF', CLU_CUTOFF, ...
                'PRECISION', PRECISION, ...
                'outDirRadius', outDirRadius, ...
                'outDirAgents', outDirAgents, ...
                'outDirClusters', outDirClusters ...      
        );
    end
  
    clusters_onefile(params);    
    truthradius_onefile(params);
    agents_onefile(params);
    
end

toc;