%% Conf parameters new Model with Anna

%%%%%%%%%%%%%

% GLOBAL Conf

simName = 'alpha-R-k-A-B-tau-sigma';
dumpDir = 'dump/';

VIDEO = 0;
DEBUG = 0;
DUMP = 1;
COMPUTATION = 2; % 0-local, 1-parallel, 2-LSF


%%%%%%%%%%%%%

% MODEL Conf

nRuns = 10;             % Number of simulation runs with same param set

dts = [0.01];           % time_step
t_ends = [30];          % running time

n_agents = [100];       % number of agents

ideas_space_sizes = [1];% size of ideas space
ideas_space_dims = [2]; % dimension of ideas space

% If A = B repulsion and attraction nullify
% If A > B and there is some initial speed -> speedy bounce is possible at boundaries

% alpha is between 0 and 1. 
%
% 0 -> velocity is just from the average within the radius;
% 1-> just your own

% ks the bigger the less groups

% VELOCITY 
alphas = [0.1:0.1:1];       	% weighting of velocity terms
Rs     = [0.1:0.1:1];       	% cut-off radius

% ATTRACTIVE AND REPULSIVE FORCES

ks     = [0.001, 0.01, 0.1, 1];           % Power of distance in force term

As     = [0:0.5:2];           % Constant in attractive force term
d0s    = [1];       	% Express the range of the interaction force (exponent divisor)

Bs     = [0:0.5:2];           % Constant in repulsive force term
d1s    = [1];       	% Express the range of the interaction force (exponent divisor)


% HOW EASY IS TO FIND THE TRUTH (
taus   = [1:2:10];     		% coupling coefficient (divisor)

% WHITE NOISE
sigmas = [0:0.1:1];       	% Std. deviation of white noise term

% INITIIAL VELOCITIES OF SCIENTISTS
vScalings = [1];     	% Scaling factor for initial (random) velocities

% INITIAL POSITIONS OF SCIENTISTS
nClusters = [0];    	% number of clusters of the initial positions
clusterTightness = [0.25]; % Tightness of clusters

% TRUTH POSITION
% Generate Truth Vector for 2D Truth

hGrid = [0:0.25:1];
vGrid = [0,0,0,0,0];
for i=2:numel(hGrid)
    vGrid = [ vGrid repmat(hGrid(i),1,5)];
end
truths = [repmat(hGrid,1,5); vGrid];

%truths = [0.5;0.5];

% BOUNDARY CONDITIONS
bBounce = 0;
bStop   = 1;
bTorus  = 2;

boundaryCondition = bBounce;

save(simName);

