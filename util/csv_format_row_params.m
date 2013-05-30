function [ row ] = csv_format_row_params( name, filename, run, simTime, params, truth )
    
    row = [ '"' name '","' filename '","'  num2str(run) '","' simTime '"' ] ;

    % There is no matching between headers and data if I use the for loop
    %cellParams = struct2cell(params);
    %for i = 1:length(cellParams)
    %    row = [row, ',' num2str(cellParams{i})];
    %end
    %row = [ row ',' num2str(truth(1,1)) ',' num2str(truth(2,1)) ',' num2str(finalConv)];
    
    row = [ row ',' num2str(params.dt) ',' num2str(params.t_end) ',' num2str(params.n_agents)];
    
    row = [ row ',' num2str(params.ideas_space_size) ',' num2str(params.ideas_space_dim) ',' num2str(params.alpha)];
    
    row = [ row ',' num2str(params.R) ',' num2str(params.k) ',' num2str(params.A)];
    
    row = [ row ',' num2str(params.d0) ',' num2str(params.B) ',' num2str(params.d1)];
    
    row = [ row ',' num2str(params.tau) ',' num2str(params.sigma) ',' num2str(params.v_scaling)];
    
    row = [ row ',' num2str(params.nof_cluster) ',' num2str(params.clusterTightness)];
    
    row = [ row ',' num2str(truth(1,1)) ',' num2str(truth(2,1))];
end

