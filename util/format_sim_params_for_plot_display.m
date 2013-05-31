function [row] = format_sim_params_for_plot_display(simname, simfile, params) 


    row = [ 'sim: ' simname ];
    
    row = [ row ' file: ' simfile];

    % row = [ ' dt: ' num2str(params.dt) ', steps: ' num2str(params.steps) ', agents: ' num2str(params.nAgents)];
    
    % row = [ row ', size: ' num2str(params.iss) ];
    
    row = [ row ' R: ' num2str(params.R) ' alpha: ' num2str(params.alpha) ];
    
    row = [ row ' tau: ' num2str(params.tau) ' sigma: ' num2str(params.sigma) ' v0: ' num2str(params.v_scaling)];
    
    % row = [ row ', A:' num2str(params.A) ', B: ' num2str(params.B) ];
    
    % row = [ row ', d0: ' num2str(params.d0) ', d1: ' num2str(params.d1)];
    
    % row = [ row ', c0: ' num2str(params.nClusters) ',ct0: ' num2str(params.clusterTightness) ', k: ' num2str(params.k)];
    
    row = [ row ' truth: ' num2str(params.truth(1,1)) ';' num2str(params.truth(2,1))];
end