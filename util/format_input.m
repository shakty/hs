function [ out ] = format_input( in )

if (size(in,1) > 1) %% that is the truth param
    out = '';
    for i=1:size(in,2)
        out = [ out '[' num2str(in(1,i)) ';' num2str(in(2,i)) '];' ];
    end
    %sizeIn = size(in);
    %out = ['[' num2str(sizeIn(1,1)) 'x' num2str(sizeIn(1,2)) '] ' num2str(in([2 1])) '...' num2str(in([end-1, end]))];
else

    in_datatype = class(in);

    switch in_datatype
        case 'char'
            out = in;
        case 'double'
            out = num2str(in);
        otherwise
            out = in;
            %error(['Unknown type ' in_datatype]);
    end
end