function [ row ] = csv_format_row_params( name, nameidx, run, t, id, pos)
    
    row = [ '"' name '", "' nameidx '", "' num2str(run) '", ' num2str(t) ] ;
    row = [row, ', ' num2str(id), ', ' num2str(pos(1,1)) ', ' num2str(pos(2,1))];
 
end

