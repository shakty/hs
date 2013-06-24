function [ output_args ] = compute_avg_stats( global_count_sum, global_count_sumsquared, CI_INT )

    nFiles = length(global_count_sum(1));
    df = nFiles - 1;

    % Computing global stats  
    t_count_avg = global_count_sum / nFiles; 
    t_count_sd = sqrt(((global_count_sumsquared - ((global_count_sum).^2 / N))) / df);
    t_count_se = t_count_sd / sqrt(nFiles);  
    t_count_ci = t_count_se * tquant(CI_INT, nFiles-1);
    
    t_cover_avg = global_coverage_sum / nFiles; 
    t_cover_sd = sqrt(((global_coverage_sumsquared - ((global_coverage_sum).^2 / N))) / df);
    t_cover_se = t_cover_sd / sqrt(nFiles);  
    t_cover_ci = t_cover_se * tquant(CI_INT, nFiles-1);
    
    t_cover_cum_avg = global_coverage_cum_sum / nFiles; 
    t_cover_cum_sd = sqrt(((global_coverage_cum_sumsquared - ((global_coverage_cum_sum).^2 / N))) / df);
    t_cover_cum_se = t_cover_cum_sd / sqrt(nFiles);  
    t_cover_cum_ci = t_cover_cum_se * tquant(CI_INT, nFiles-1);
    
    t_speed_avg = global_speed_sum / nFiles; 
    t_speed_sd = sqrt(((global_speed_sumsquared - ((global_speed_sum).^2 / N))) / df);
    t_speed_se = t_speed_sd / sqrt(nFiles);  
    t_speed_ci = t_speed_se * tquant(CI_INT, nFiles-1);
    
    t_move_avg = global_move_sum / nFiles; 
    t_move_sd = sqrt(((global_move_sumsquared - ((global_move_sum).^2 / N))) / df);
    t_move_se = t_move_sd / sqrt(nFiles);  
    t_move_ci = t_move_se * tquant(CI_INT, nFiles-1);
    
    t_size_avg = global_size_sum / nFiles; 
    t_size_sd = sqrt(((global_size_sumsquared - ((global_size_sum).^2 / N))) / df);
    t_size_se = t_size_sd / sqrt(nFiles);  
    t_size_ci = t_size_se * tquant(CI_INT, nFiles-1);

    t_fromtruth_avg = global_fromtruth_sum / nFiles; 
    t_fromtruth_sd = sqrt(((global_fromtruth_sumsquared - ((global_fromtruth_sum).^2 / N))) / df);
    t_fromtruth_se = t_fromtruth_sd / sqrt(nFiles);  
    t_fromtruth_ci = t_fromtruth_se * tquant(CI_INT, nFiles-1);
end

