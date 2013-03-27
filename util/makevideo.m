function makevideo( fileIn, fileOut, MPEG )

    DEBUG = 0;

    colors = {'magenta','yellow','black', 'cyan', 'red', 'green', 'blue'};
    
    load(fileIn);
   
    allStepsAgents = dump.agents;
    truth = dump.truth;
    
    % Creating a string with the description of the parameters
    paramString = ['File: ' fileIn];
    paramString = [paramString create_params_string(dump.parameters, dump.truth)];
    

    %% Video Plotting

    
    xlim([0 1]);
    ylim([0 1]);
    

    % Get the handle of the figure
    %h = figure();

    if (MPEG)
        % Prepare the new file.
        vidObj = VideoWriter(fileOut);
        open(vidObj);
    end

    %A = allStepsAgents;  
    %figure
    close all
    for j=1:size(allStepsAgents,3)

        if (j < 4)
            continue;
        end

        agents = allStepsAgents(:,:,j);

        % PLOT red crosses
        plot(agents(1,:),agents(2,:),'rx');     

        % PLOT COLORED NUMBERS
        %points = arrayfun(@(x) {[ '\color{' colors{mod(x,length(colors))+1} '}' int2str(x)]}, 1:length(agents));
        %text(agents(1,:),agents(2,:), points');

        % PLOT BLACK NUMBERS
        %text(agents(1,:),agents(2,:), num2str([1:length(agents)]'));

        hold on;
        %plot(agents_average(1),agents_average(2),'bo');
        plot(truth(1),truth(2),'go');
        title(['T: ' int2str(j) ' ' paramString]);

        hold off;

        % LIMITS for plotting red crosses
        xlim([0 dump.parameters.iss]);
        ylim([0 dump.parameters.iss]);

        % LIMITS when plotting numbers
        %axis([0 dump.parameters.iss 0 dump.parameters.iss])

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

        %if (j<500)
            pause(0.01);
        %else
        %    pause(0.2);
        %end
        % no need for pause when plotting numbers

        %if (j ~= size(allStepsAgents,3)-1)
        %    clf
        %end
    end      
 
%%    
    if (MPEG)
        % Close the video file.
        close(vidObj)
    end
end

