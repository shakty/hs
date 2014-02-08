function [agents_pos] = initial_pos_clustered(nof_clusters, ...
    clusterTightness, n_agents, ideas_space_size, ideas_space_dim, ...
    clustersInCircleOfRadius, truth)

    % nof_clusters can be:
    %
    %   -A (i) a number, or (ii) 1D array -> indicates the number of clusters
    %   -B (i) a 2D, or (ii) 3D array -> indicates the positions of the centers

    % If A, clustersInCircleOfRadius is used. It can be:
    %
    %  -1 -> clusters randomly placed
    %  any value > 0 -> all clusters will be placed on the radius
    
    
    % To be backward compatible we need to do a bit more checking.
    sizeNC = size(nof_clusters);
    
    % centers is just a number instead of an array of centers
    if (sizeNC(1,1) == 1 && nof_clusters(1,1) == 0)
        % Agents are positioned randomly across the whole idea space
        agents_pos = ideas_space_size.*rand(ideas_space_dim,n_agents);    
    
    else
        
        if (sizeNC(1,1) == 1)
            % N centers are created.
            
            if (clustersInCircleOfRadius == -1)
                
                % Randomly placed.
                centers = ideas_space_size .* rand(ideas_space_dim, nof_clusters(1,1));
            else
                
                % On a radius from truth.
                centers = circle(truth(1,1), truth(2,1), ... 
                    clustersInCircleOfRadius, nof_clusters(1,1));
            end
            
        else
            % N centers were already passed.
            centers = nof_clusters;
        end

        nC = size(centers, 2);
        
        % Uniformly distributed random radii of the clusters.
        radii = (ideas_space_dim / nC) .* rand(1, nC);

        agents_in_clusters = zeros(1, nC);

        % "assign" agents to a random cluster
        for i = 1 : n_agents
            cIdx = mod(i, nC);
            if (cIdx == 0) 
                cIdx = nC;
            end
           agents_in_clusters(cIdx) = agents_in_clusters(cIdx) + 1; 
        end

        agents_pos = [];

        for i = 1 : nC
            % for every cluster, randomly (normally distributed) distribute the
            % corresponding agents around the center of the cluster, more or
            % less within the calculated radius
            cluster_pos = clusterTightness * radii(i) * ...
                            randn(ideas_space_dim, agents_in_clusters(i)) ... 
                            + repmat(centers(:,i), 1, agents_in_clusters(i));
            
            agents_pos = [agents_pos, cluster_pos];                   
        end
        
        % in case agents are positioned outside of feasible area: 
        % position them on the border
        dimHa = agents_pos(1,:) < 0;
        dimHb = agents_pos(1,:) > ideas_space_size;
        dimH = bitor(dimHa, dimHb);
        agents_pos(1,dimH) = round(agents_pos(1,dimH));
        dimVa = agents_pos(2,:) > ideas_space_size;
        dimVb = agents_pos(2,:) < 0;
        dimV = bitor(dimVa, dimVb);
        agents_pos(2,dimV) = round(agents_pos(2,dimV));
    end
end


function h = circle(x, y, radius, howmany)  
    th = 0 : 2*pi / howmany : 2*pi;
    xunit = radius * cos(th) + x;
    yunit = radius * sin(th) + y;
    h = [xunit ; yunit];
end