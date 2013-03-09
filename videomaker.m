%% Perform Cluster Analysis on the dump

close all
clc
clear

path(path,'../util/');

VIDEO = 1;
MPEG  = 1;
DEBUG = 0;

DUMPDIR = 'dump/';
VIDEODIR = 'videos/';

%'the_loop-2013-3-8-16-13/';
%MYDIR = 'circle_maybe-2013-3-8-13-22/';
%MYDIR = 'many_little_clusters_do_not_find_truth_mod_a_little_bigger_clusters-2013-3-8-17-29/'
%MYDIR = 'few_big_groups_do_not_find_truth_in_a_smaller_space_cluster_together-2013-3-8-18-28/';

%MYDIR = 'few_big_groups-DIM-vs-ALPHA/';

% WITH THE BUG (full) 'the_loop-2013-3-8-16-13/'
% WITH THE BUG (withtout A and B) 'the_loop-2013-3-9-16-23/'
% WITHOUT THE BUG (no alpha - R): 'the_loop-2013-3-9-16-29/'
% the_loop-2013-3-9-16-42-a/


MYDIR = 'R-alpha-noA-noB-2013-3-6-20-8/'; 

MYDIR = 'tests/the_loop-2013-3-9-22-40/'

dumpDir = [DUMPDIR MYDIR]; 

%dumpDir = [DUMPDIR '2011-5-30-15-48/'];
%dumpDir = [DUMPDIR '2011-5-30-15-54/'];
%dumpDir = [DUMPDIR '2011-5-30-15-56/'];
%dumpDir = [DUMPDIR '2011-5-30-16-3/'];

files = dir(dumpDir);
fileIndex = find(~[files.isdir]);

%%

paramOrder = {' agents=',' A=',' B=',' k=',' d0=','d1=',...
    ' alpha=',' tau=',' R=',' sigma=',' velocity scaling=',...
    ' n. of clusters='};

display(['Files Found: ' num2str(length(fileIndex))]);

if (~length(fileIndex) > 1)
    display('Empty dir.');
end

%%

colors = {'magenta','yellow','black', 'cyan', 'red', 'green', 'blue'};
%%

%for i = 1:length(fileIndex)
    
    i=1
    
    append = files(fileIndex(i)).name;
    fileName = [dumpDir, append];
    
     % We load only .mat
     [PATH,NAME,EXT] = fileparts(fileName);
     if (~strcmpi(EXT,'.mat')) 
        continue;
     end
    
    load(fileName);
    
    simParameters = cell2mat(struct2cell(dump.parameters));
    allStepsAgents = dump.agents;

    % Creating a string with the description of the parameters
    paramString = ['File: ' append];
    paramString = [paramString create_params_string(dump.parameters, dump.truth)];
    
    %for j=1:size(paramOrder,2)
     %   paramString = [paramString ';' paramOrder{j} ...
      %      int2str(simParameters(j)) ];
    %end
    
    paramString
    
    %% Video Plotting

    
    xlim([0 1]);
    ylim([0 1]);
    
         
    
    if (VIDEO)
        
        % Get the handle of the figure
        %h = figure();

        if (MPEG)
            % Prepare the new file.
            vidObj = VideoWriter([VIDEODIR append '.avi']);
            open(vidObj);
        end
        
        %A = allStepsAgents;  
        %figure
        close all
        for j=1:size(allStepsAgents,3)
            truth = dump.truth;
            agents = allStepsAgents(:,:,j);
            
            % PLOT red crosses
            %plot(agents(1,:),agents(2,:),'rx'); 
            
            
            % PLOT COLORED NUMBERS
            %points = arrayfun(@(x) {[ '\color{' colors{mod(x,length(colors))+1} '}' int2str(x)]}, 1:length(agents));
            %text(agents(1,:),agents(2,:), points');
            
            % PLOT BLACK NUMBERS
            text(agents(1,:),agents(2,:), num2str([1:length(agents)]'));
            
            hold on;
            %plot(agents_average(1),agents_average(2),'bo');
            plot(truth(1),truth(2),'go');
            title(['T: ' int2str(j) ' ' paramString]);
            
            hold off;
            
            % LIMITS for plotting red crosses
            %xlim([0 dump.parameters.iss]);
            %ylim([0 dump.parameters.iss]);
            
            % LIMITS when plotting numbers
            axis([0 dump.parameters.iss 0 dump.parameters.iss])
            
            if (MPEG)
                % Get the very last frame
                currFrame = getframe();
                writeVideo(vidObj,currFrame);
                %A(j)=getframe();
            end 

            if (DEBUG)
                % Debugging: Calculate energy and control...
                % whether it is constant
                E=0;
                for ii=1:n_agents
                    E=E+0.5*(v(:,ii))'*v(:,ii);
                end
                energy=num2str(E);
                legend(energy);
            end

            % no need for pause when plotting numbers
            %pause(0.01);
            if (j ~= size(allStepsAgents,3)-1)
                clf
            end
            
        end
    end
 
%%    
    % sscanf(paramString,'%s')
    if (MPEG)
        % Close the video file.
        close(vidObj)
        
        %movie2avi(A, ['videos/' append '.avi'],'fps',60);
        %save movie1.mat A
    end
    
%end    
