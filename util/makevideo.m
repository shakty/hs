function makevideo( fileIn, MPEG, fileOut, plottype, SHOW_POTENTIAL)
    close all;
    
    % PLOT TYPE
    plot_cross = 0;
    plot_number = 1;
    plot_number_color = 2;
    plot_arrow = 3;
    % Attraction to truth force type
    attr_zero = 0;
    attr_const = 1;
    attr_linear = 2;
    attr_expo = 3;
    attr_millean_arena = 4;
    attr_hard_to_find = 5;
    attr_wide_funnel = 6;
    attr_gentle_landing = 7;

    if nargin < 4
        plottype = plot_cross;
    
    elseif nargin < 5
        SHOW_POTENTIAL = 0;
    end
    
    DEBUG = 0; % not used for now
 
    colors = {'magenta','yellow','black', 'cyan', 'red', 'green', 'blue'};
    
    load(fileIn);
    
    
    
   
    %allStepsAgents = dump.agents;
    agents = dump.agents;
    v = dump.agentsv;
    truth = dump.truth;
    
    % not used for now
    % Creating a string with the description of the parameters
    %paramString = ['File: ' fileIn];
    %paramString = [paramString create_params_string(dump.parameters, dump.truth)];
    

    %% Video Plotting

    % Showing the potential of the attraction to truth
    if (SHOW_POTENTIAL)
        attrtype = dumps.parameters.attrtype;
    
        switch (attrtype)

            case attr_zero
            % NO TRUTH
            ths = @(x) (zeros(2, n_agents));

            case attr_const
            % TRUTH Constant
            ths = @(x) (repmat(truth,1,length(x))-x)./tau.*(repmat(tau./colnorm(repmat(truth,1,length(x))-x,2),2,1));

            case attr_linear
            % TRUTH Linear
            ths = @(x) (repmat(truth,1,length(x))-x)./tau;

            case attr_expo
            % TRUTH Decaying exponentially (EXP)
            SIGMA = 1;
            ths = @(x) (repmat(truth,1,length(x))-x)./tau.*repmat(exppdf(colnorm(repmat(truth,1,length(x)) - abs((repmat(truth,1,length(x))-x)),2),SIGMA),2,1);   

            case attr_millean_arena
            % Millean Arena (NORMAL)
            POS = 3; SIGMA = 0.05;
            ths = @(x) (repmat(truth,1,length(x))-x)./tau.*repmat(normpdf(colnorm(repmat(truth,1,length(x))-x,2),(DIAG -norm(truth))/POS,SIGMA),2,1);

            case attr_hard_to_find
            % Hard to Find (NORMAL)
            POS = 100; SIGMA = 0.02;
            ths = @(x) (repmat(truth,1,length(x))-x)./tau.*repmat(normpdf(colnorm(repmat(truth,1,length(x))-x,2),(DIAG -norm(truth))/POS,SIGMA),2,1);

            case attr_wide_funnel
            % Wide Funnel to Truth (LOG-NORMAL)
            POS = 1; SIGMA = 3;
            ths = @(x) (repmat(truth,1,length(x))-x)./tau.*repmat(lognpdf(colnorm(repmat(truth,1,length(x))-x,2),(DIAG -norm(truth))/POS,SIGMA),2,1);

            case attr_gentle_landing
            % Gentle Landing to truth
            POS = 0; SIGMA = 0.2;
            ths = @(x) (repmat(truth,1,length(x))-x)./tau.*repmat(normpdf(colnorm(repmat(truth,1,length(x))-x,2),POS,SIGMA),2,1);   

        end
        
        if (attrtype > 1)
            PRECISION = 0.01;
            a = [0:PRECISION:1;0:PRECISION:1];
            [X,Y] = meshgrid(a(1,:), a(1,:));
            Z = zeros(size(X));
            for i=1:length(X)
                potential_grid = [X(i,:) ; Y(i,:)];
                forces = ths(potential_grid);
                Z(i,:) = colnorm(forces,2);
            end
            contour(X,Y,Z);
        end
    end
    
    
    xlim([0 1]);
    ylim([0 1]);
    

    % Get the handle of the figure
    %h = figure();

    if (MPEG)
        % Prepare the new file.
        vidObj = VideoWriter(fileOut);
        open(vidObj);
    end

    for j=1:size(allStepsAgents,3)
   
        switch (plottype)
        
            case plot_cross
            % PLOT red crosses
            plot(agents(1,:),agents(2,:),'rx');     
        
            case plot_number
            % PLOT BLACK NUMBERS
            text(agents(1,:),agents(2,:), num2str([1:length(agents)]'));
            plot(exp(agents(1,1)));
        
            case plot_number_color
            % PLOT COLORED NUMBERS
            points = arrayfun(@(x) {[ '\color{' colors{mod(x,length(colors))+1} '}' int2str(x)]}, 1:length(agents));
            text(agents(1,:),agents(2,:), points');
            
            case plot_arrow
            % PLOT VELOCITY ARROWS
            quiver(agents(1,:),agents(2,:),v(1,:),v(2,:));
            

        end
        hold on;
        
        if (SHOW_POTENTIAL && attrtype > 1)
            [C, h] = contour(X,Y,Z);
            alpha(.5);
        end
        
        plot(truth(1),truth(2),'go');
        
        hold off;
            
        xlim([0 ideas_space_size]);
        ylim([0 ideas_space_size]);
        
        if (MPEG)
            % Get the very last frame
            currFrame = getframe();
            writeVideo(vidObj,currFrame);
        end 
        
        if (DEBUG)
            %Debugging: Calculate energy and control whether it is constant
            E=0;
            for i=1:n_agents
                E=E+0.5*(v(:,i))'*v(:,i);
            end
            energy=num2str(E);
            legend(energy);
        end
        
        % no need for pause for complicated plots
        if (params.plottype == plot_cross)
            pause(0.01);
        end
       
    end      
 
%%    
    if (MPEG)
        % Close the video file.
        close(vidObj)
    end
end

