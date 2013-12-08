function write_csv_headers( filename, headers )

%% write the header string to the file

%turn the headers into a single comma seperated string if it is a cell
%array, 
header_string = ['"' headers{1} '"'];
for i = 2:length(headers)
    header_string = [header_string,', "', headers{i} '"'];
end

%write the string to a file
fid = fopen(filename,'w');
fprintf(fid,'%s\n', header_string);
fclose(fid);

end

