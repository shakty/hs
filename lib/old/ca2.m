%% Perform Cluster Analysis on the dump

close all
clc
clear

path(path,'analysis');
path(path,'util');


LASTFRAME   = 0;
VIDEO       = 1;
DEBUG       = 0;

DUMPDIR = 'dump/';

% Do Not forget trailing slashes /
dumpDir = [DUMPDIR '30A_Taus_bouncing-2011-6-26-13-17/']; 

files = dir(dumpDir);
fileIndex = find(~[files.isdir]);

if (numel(fileIndex)==0)
    error('No files in the directory'); 
end

nCells = 256;

grid = zeros(256);

for i = 1:length(fileIndex)
    append = files(fileIndex(i)).name;
    fileName = [dumpDir, append];
    load(fileName); % load a struct variable called dump
    
    simName = dump.name;
    simSim = dump.sim;
    simRun = dump.run;
    
    simParameters = dump.parameters;
    allStepsAgents = dump.agents;
    
    paramNames = fieldnames(simParameters);
    
    % Creating a string with the description of the parameters
    paramString = ['File: ' append ' '];
    for j=1:numel(paramNames)
        name = paramNames{j};
        value = simParameters.(paramNames{j});
        paramString = [paramString name ' ' int2str(value) '; '];
        %paramString = [paramString int2str(j)];
    end
    
    %display([simName ': ' int2str(simSim) '.' int2str(simRun)]);
    %display(paramString);
    %paramString
    grid = grid+countAgents(allStepsAgents(:,:,end),256);
end

%save([dumpDir 'grid.mat'],'grid');
%createAndSaveFigure(grid,[dumpDir 'grid.png']);


%imagesc(grid);  

if (LASTFRAME)
    figure
    lastFrame = allStepsAgents(:,:,end);
    plot(lastFrame(1,:),lastFrame(2,:),'rx');
    title(['T: ' int2str(j) ' ' paramString]);

    hold off;
    xlim([0 1]);
    ylim([0 1]);
end

if (VIDEO)
    figure
    for j=1:size(allStepsAgents,3)
        truth = dump.truth;
        agents = allStepsAgents(:,:,j);
        plot(agents(1,:),agents(2,:),'rx'); 
        hold on;
        %plot(agents_average(1),agents_average(2),'bo');
        plot(truth(1),truth(2),'go');
        title(['T: ' int2str(j) ' ' paramString]);

        hold off;
        xlim([0 1]);
        ylim([0 1]);

        pause(0.01);
    end
 end
 