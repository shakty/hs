%% Perform Cluster Analysis on the dump

close all
clc
clear

%path(path,'../../FUZZCLUST');

VIDEO = 1;
MPEG  = 1;
DEBUG = 0;

DUMPDIR = 'dump/';

%'the_loop-2013-3-8-16-13/';
%MYDIR = 'circle_maybe-2013-3-8-13-22/';
MYDIR = 'many_little_clusters_do_not_find_truth_mod_a_little_bigger_clusters-2013-3-8-17-29/'

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


%%

for i = 1:length(fileIndex)
    
    
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
            vidObj = VideoWriter(['videos/' append '.avi']);
            open(vidObj);
        end
        
        %A = allStepsAgents;  
        %figure
        close all
        for j=1:size(allStepsAgents,3)
            truth = dump.truth;
            agents = allStepsAgents(:,:,j);
            plot(agents(1,:),agents(2,:),'rx'); 
            hold on;
            %plot(agents_average(1),agents_average(2),'bo');
            plot(truth(1),truth(2),'go');
            title(['T: ' int2str(j) ' ' paramString]);
            
            hold off;
            xlim([0 dump.parameters.iss]);
            ylim([0 dump.parameters.iss]);

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

            pause(0.01);
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
    
end    
