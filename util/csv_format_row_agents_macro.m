function [ row ] = csv_format_row_agents_macro(obj, simname, t)
    
    row = sprintf('"%s",%u,%u,%u,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f', ... 
        simname, ...
        obj.simnameidx, ...
        obj.run, ... 
        t, ...       
        obj.avgcoverage, ...
        obj.cumcoverage, ...
        obj.mean_agents_speed, ...
        obj.sd_agents_speed, ...        
        obj.mean_agents_movs, ...
        obj.sd_agents_movs, ...
        obj.mean_agents_fromtruth, ...
        obj.sd_agents_fromtruth, ...
        obj.pairwise_dist_mean, ...
        obj.pairwise_dist_sd ...        
    );
      
end