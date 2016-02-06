function obj=isp_unpickle(vec)
%ISP_UNPICKLE  Recover original Matlab object from pickled column vector.
%
% SYNTAX
%     obj = isp_unpickle(vec)
%
% DESCRIPTION
%   Using isp_pickle, a Matlab object build from the basic Matlab types
%   can be converted to a vector of doubles for easy storage in e.g. a
%   database. isp_unpickle recovers the original object from the vector.
%
% INPUT
%   vec:
%     A column vector of doubles as returned by isp_pickle.
%
% OUTPUT
%   obj:
%     A copy of the object passed to isp_pickle.
%     logical, char, the integer types, cells and structs). Sparse
%     matrices will be converted to full matrices.
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


    % Ensure column vector
    vec = vec(:);

    datatypes = {'double', 'single', 'logical', 'char', 'int8', 'uint8', ...
                 'int16', 'uint16', 'int32', 'uint32', 'int64', 'uint64', ...
                 'cell', 'struct'};

    datafunc = {@double, @single, @logical, @char, @int8, @uint8, ...
                 @int16, @uint16, @int32, @uint32, @int64, @uint64};

    obj=isp_unpickle_helper(vec, 1);


    function [obj, dataidx]=isp_unpickle_helper(vec, dataidx)

        typeNo = vec(dataidx);
        nDim = vec(dataidx+1);
        dim=vec(dataidx+2:dataidx+nDim+1)';
        dataidx=dataidx + nDim+2;
        if typeNo <= 12
            dataend = dataidx+prod(dim)-1;
            obj=reshape(datafunc{typeNo} (vec(dataidx:dataend)), dim);
            dataidx = dataend+1;
        elseif typeNo == 13 % cell
            obj = cell(dim);
            for n=1:numel(obj)
                [obj{n}, dataidx] = isp_unpickle_helper(vec, dataidx);
            end
        elseif typeNo == 14 % struct
            obj = struct([]);
            nFields = vec(dataidx);
            dataidx = dataidx+1;
            for m=1:nFields
                fieldnamelen = vec(dataidx);
                fieldname = char(vec(dataidx+1:dataidx+fieldnamelen)');
                dataidx = dataidx + fieldnamelen+1;
                for n=1:prod(dim)
                    [obj(n).(fieldname), dataidx] = isp_unpickle_helper(vec, dataidx);
                end
            end
            obj=reshape(obj, dim);
        else
            error('This is not supposed to happen ...')
        end

    end

end
