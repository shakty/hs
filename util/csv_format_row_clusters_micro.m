function [ row ] = csv_format_row_clusters_micro( obj)
    
    for i=1:length(obj.clusters_speed)
        row = [ '"' obj.name '","' obj.nameidx '","'  num2str(obj.run) '",' num2str(obj.t) ] ;
        row = [row, '","' num2str(obj.size(i)), '","' num2str(obj.speed(i)) ];
        row = [row, '","' num2str(obj.move(i)), '","' num2str(obj.fromtruth(i)) ];
    end
    
end

