%% Saves simulations into properly formatted CSV files

close all;
clear;
clc;

%% Add other directories to path
path(path,'util/'); % Help functions

%% Param

headers_data = {
    'sim', ...
    'run', ...
    't', ...
    'id', ...
    'x', ...
    'y'};

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


DUMPDIR = 'dump/';


%% Load Data


simName = 'few_big_groups-DIM-vs-ALPHA';


dumpDir = [DUMPDIR simName '/'];


files = dir(dumpDir);
fileIndex = find(~[files.isdir]);

if (isempty(fileIndex))
    error('Invalid Directory Selected');
end

%% Date and Time

mytimestamp = datestr ( datevec ( now ), 0 );


%% Write headers

paramFileName = [dumpDir 'params.csv'];
dataFileName = [dumpDir 'data.csv'];

write_csv_headers(paramFileName, headers_params);
write_csv_headers(dataFileName, headers_data);

%% Write content into csv files

fidParam = fopen(paramFileName,'a');
fidData = fopen(dataFileName,'a');

% Load all parameters matrices in one
for i = 1:length(fileIndex)
    
    append = files(fileIndex(i)).name;
    fileName = [dumpDir, append];
    load(fileName);
    
    param_string = csv_format_row_params(simName, dump.run, mytimestamp, dump.parameters, dump.truth, dump.conv);
    
    
    % append the param string to the file
    fprintf(fidParam,'%s\n', param_string);

    t = 1;
    for z=1:size(dump.agents, 3)
       
        agentpos = dump.agents(:,:,z);
        for id = 1:size(agentpos,2)     
            data_string = csv_format_row_data(simName, dump.run, t, id, agentpos(:,id));
            fprintf(fidData,'%s\n', data_string);
        end
        t = t + 1;
    end
    
    
end


fclose(fidParam);
fclose(fidData);


%%


X = agentpos';
Y = pdist(X,'euclidean');
Z = linkage(Y,'average');
%[H,T] = dendrogram(Z,'colorthreshold','default');
[H,T] = dendrogram(Z,'colorthreshold',0.1);
set(H,'LineWidth',2)
 