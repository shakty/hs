function struct2File( s, dumpDir, simName)
%Write struct to text file. The data in the struct can be both numbers and
%strings. 

    fileName = [ dumpDir simName '/params.txt'];

    fields = fieldnames(s)';

    fid = fopen(fileName,'wt');
    if (fid==-1)
        error(['Could not open fileName ' fileName]);
    end


    fprintf(fid, '%s\n',simName);
    fprintf(fid, '------------------------------------\n');

    for k=1:numel(fields)
        field = fields{k};
        value = s.(field);
        value = format_input(value);
        fprintf(fid, '%s = %s\n', field, value);
    end

    fclose(fid);
    
   
end