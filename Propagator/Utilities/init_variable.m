function [ var ] = init_variable( nrows, mcols, zslices, datatype, flag )
%[ var ] = make_var( nrows, mcols, datatype )
%   Function to initialize a variable of size [nrows,mcols,zslices] of
%   data type datatype. Data Types supported are single, double, and uint8.
%   If flag is set to 1, initialize to ones matrix, else make it zeros.

if nargin < 5
    flag = 0;
end

switch datatype
    case 'single'
        if flag
            var = single(ones(nrows,mcols,zslices));
        else
            var = single(zeros(nrows,mcols,zslices));
        end
        
    case 'double'
        if flag
            var = double(ones(nrows,mcols,zslices));
        else
            var = double(zeros(nrows,mcols,zslices));
        end
    case 'uint8'
        if flag
            var = uint8(ones(nrows,mcols,zslices));
        else
            var = uint8(zeros(nrows,mcols,zslices));
        end
    otherwise
        error('Data Types Supported are single, double, and uint8');
end




end

