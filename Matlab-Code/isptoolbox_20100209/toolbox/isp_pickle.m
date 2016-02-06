function vec=isp_pickle(obj)
%ISP_PICKLE  Convert Matlab data to column vector
%
% SYNTAX
%     vec = isp_pickle(obj)
%
% DESCRIPTION
%   Converts Matlab data to a vector of doubles for storage in e.g. a
%   database.
%
% INPUT
%   obj:
%     A Matlab object build from simple Matlab types (double, single,
%     logical, char, the integer types, cells and structs). Sparse
%     matrices will be converted to full matrices.
%
% OUTPUT
%   vec:
%     A column vector of doubles from which isp_unpickle can recover the
%     original object, 'obj'.
%
% EXAMPLE
%   Create a struct, and 'pickle' and 'unpickle' it:
%     >> a1=struct('f', 34, 'q', 'fgdsa', 'h', {{34,1.4, 8}})
%     a1 =
%         f: 34
%         q: 'fgdsa'
%         h: {[34]  [1.4000]  [8]}
%     >> b=isp_pickle(a)
%     b =
%        14.0000
%         2.0000
%         1.0000
%         ...
%         8.0000
%     
%     >> a2=isp_unpickle(b)
%     a2 =
%         f: 34
%         q: 'fgdsa'
%         h: {[34]  [1.4000]  [8]}
%
% SEE ALSO
%   isp_unpickle.
%
% HISTORY
%   Created by Jesper H. Jensen January 2008.

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2 as published
% by the Free Software Foundation.


    datatypes = {'double', 'single', 'logical', 'char', 'int8', 'uint8', ...
                 'int16', 'uint16', 'int32', 'uint32', 'int64', 'uint64', ...
                 'cell', 'struct'};

    typeNo = find(strcmp(class(obj), datatypes));

    vec = [];

    if isempty(typeNo)
        error(['Unknown data type ' class(obj) '.'])
    elseif typeNo <= 12
        %vec(1:end+ndims(obj)+2+numel(obj), 1) = [typeNo; ndims(obj); size(obj)'; obj(:)];
        vec = double([typeNo; ndims(obj); size(obj)'; obj(:)]);
    elseif typeNo == 13 % cell
        tmp = cell(numel(obj),1);
        for n=1:numel(obj)
            tmp{n} = isp_pickle(obj{n});
        end
        vec = [typeNo; ndims(obj); size(obj)'; cat(1, tmp{:})];

    elseif typeNo == 14 % struct
        fina = fieldnames(obj)';
        vec = [typeNo; ndims(obj); size(obj)'; length(fina)];
        tmp = cell(numel(obj),1);
        for fn=fina
            for n=1:numel(obj)
                tmp{n} = isp_pickle(obj(n).(fn{1}));
            end
            vec = [vec; length(fn{1}); double(fn{1})'; cat(1, tmp{:})];
        end
    else
        error('This is not supposed to happen ...')
    end

end
