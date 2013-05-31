function [ row ] = csv_format_row_clusters_micro(stepData, simname, t)
    
    for i=1:length(stepData.clusters_speed)
        row = [ '"' simname '","' stepData.simnameidx '","'  num2str(stepData.run) '",' num2str(t) ] ;
        row = [row, '","' num2str(stepData.clusters_size(i)), '","' num2str(stepData.clusters_speed(i)) ];
        row = [row, '","' num2str(stepData.clusters_move(i)), '","' num2str(stepData.clusters_fromtruth(i)) '"' ];
    end

   
    
end

