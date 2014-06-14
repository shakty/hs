function makevideo( dirIn, fileIn, MPEG, fileOut, plottype, SHOW_POTENTIAL, LIMITS)
    close all;
    
    path2file = [ dirIn fileIn ];
    
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
    
    elseif nargin < 6
        LIMITS = 1;
    end
    
    DEBUG = 0; % not used for now
 
    colnorm = @(X,P) sum(abs(X).^P,1).^(1/P);
    
    colors = {'magenta','yellow','black', 'cyan', 'red', 'green', 'blue'};
    
    load(path2file);
    
    p = dump.parameters;
    
    fprintf('\n%s\n', path2file);
    fprintf('------------------------------------\n');
    fprintf('%+15s = %2.3f\n','R',p.R);
    fprintf('%+15s = %2.3f\n','alpha',p.alpha);    
    fprintf('%+15s = %2.3f\n','tau',p.tau)
    fprintf('%+15s = %2.3f\n','sigma',p.sigma);
    fprintf('%+15s = %2.3f\n','epsilon',p.epsilon);    
    fprintf('%+15s = %2.3f\n','v_scaling',p.v_scaling);
    fprintf('%+15s = %d\n','steps',p.t_end);
    fprintf('%+15s = %d\n','nAgents',p.n_agents);
    fprintf('%+15s = %1d\n','IdeasSpace size',p.ideas_space_size);    
    fprintf('%+15s = [%2.3f:%2.3f]\n','truth',p.truth(1,1), p.truth(2,1));
    fprintf('%+15s = %d\n', 'Attr. type', p.attrtype);
    fprintf('%+15s = %d\n', 'Noise type', p.noisetype);
    fprintf('%+15s = %d\n', 'Forces on V', p.forces_on_v);
    fprintf('%+15s = %d\n','Seed', p.seed);
    fprintf('------------------------------------\n');

    allStepsAgents = dump.agents;
    truth = dump.truth;
    ideas_space_size = dump.parameters.ideas_space_size;
        
    % not used for now
    % Creating a string with the description of the parameters
    paramString = {};
    paramString{2} = create_params_string_small(dump.parameters);
    
    
    %% Video Plotting
    attrtype = dump.parameters.attrtype;
    
    switch (attrtype)

        case attr_zero
        attrName = 'Relativistic World';

        case attr_const
        attrName = 'CONST';

        case attr_linear
        attrName = 'Normal Science';
        
        case attr_expo
        attrName = 'Exp';

        case attr_millean_arena
        attrName = 'Millean Arena';
            
        case attr_hard_to_find
        attrName = 'Revolutionary Science';
        
        case attr_wide_funnel
        attrName = 'Funnel';

        case attr_gentle_landing
        attrName = 'Gentle';

    end
    
    % Showing the potential of the attraction to truth
    if (SHOW_POTENTIAL && attrtype > 1)  
        tau = dump.parameters.tau;
        colnorm = @(X,P) sum(abs(X).^P,1).^(1/P);
        DIAG = norm([ideas_space_size;ideas_space_size]);
        
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
            a = [0:PRECISION:ideas_space_size; 0:PRECISION:ideas_space_size];
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
    
    
    %xlim([0 1]);
    %ylim([0 1]);
    

    % Get the handle of the figure
    %h = figure();

    if (MPEG)
        % Prepare the new file.
        vidObj = VideoWriter(fileOut);
        open(vidObj);
    end

    for j=1:size(allStepsAgents,3)
        agents = allStepsAgents(:,:,j);
        %v = dump.agentsv(:,:,j);
        switch (plottype)
        
            case plot_cross
            % PLOT red crosses
            plot(agents(1,:), agents(2,:),'ok', ...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','r',...
                'MarkerSize',5);     
        
            case plot_number
            % PLOT BLACK NUMBERS
            cla
            text(agents(1,:),agents(2,:), num2str([1:length(agents)]'));
           
        
            case plot_number_color
            % PLOT COLORED NUMBERS
            cla
            points = arrayfun(@(x) {[ '\color{' colors{mod(x,length(colors))+1} '}' int2str(x)]}, 1:length(agents));
            text(agents(1,:),agents(2,:), points');
            
            case plot_arrow
            % PLOT VELOCITY ARROWS
            v = dump.agentsv;
            % no rescaling
            quiver(agents(1,:),agents(2,:),v(1,:,j),v(2,:,j), 0);
            %quiver(agents(1,:),agents(2,:),v(1,:,j),v(2,:,j));
            mean(abs(colnorm(v(:,:,j),2)))
        end
        
        hold on;
        
        if (SHOW_POTENTIAL && attrtype > 1)
            [C, h] = contour(X,Y,Z);
            alpha(.5);
        end        
        
        if (attrtype > 1)
            plot(truth(1),truth(2),'go');
        end
        
        hold off;
        
        if (LIMITS)
            xlim([0 ideas_space_size])
            ylim([0 ideas_space_size]);
        end
        
        % TODO move the title in the first frame
        paramString{1} = [attrName ' - T: ' int2str(j)];
        title(paramString);

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
        
        pause(0.01);
        
        % avg_v = mean(mean(abs(dump.agentsv(:,:,j))))
        if (j == 2000) 
            mov = mean(colnorm(dump.agentsv(:,:,j),2))        
            dis = mean(colnorm(agents-(repmat(truth,1, 100)),2))
        end
    end      
 
%%    
    if (MPEG)
        % Close the video file.
        close(vidObj)
    end
end

