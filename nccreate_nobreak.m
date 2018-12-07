function nccreate_nobreak(ncFile, varName, varargin))
    try
        nccreate(ncFile, varName, varargin);
    catch         
    end
end

