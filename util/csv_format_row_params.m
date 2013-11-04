function [ row ] = csv_format_row_params( name, filename, run, simTime, params, truth )
    

    row = sprintf('"%s",%u,%u,"%s",%f,%u,%u,%u,%u,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%u,%f,%f,%u,%u,%u', ... 
        name, ...
        filename, ...
        run, ... 
        simTime, ...
        params.dt, ...
        params.t_end, ...
        params.n_agents, ...
        params.ideas_space_size, ...
        params.ideas_space_dim, ...
        params.alpha, ...
        params.R, ...
        params.k, ...
        params.A, ...
        params.d0, ...
        params.B, ...
        params.d1, ...
        params.tau, ...
        params.sigma, ...
        params.v_scaling, ...
        params.nof_cluster, ...
        params.clusterTightness, ...
        truth(1,1), ...
        truth(2,1), ...
        params.noisetype, ...
        params.attrtype, ...
        params.forces_on_v ...
    );
end

