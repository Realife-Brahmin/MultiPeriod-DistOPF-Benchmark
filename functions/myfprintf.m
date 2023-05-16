function myfprintf(verbose, varargin)
    if verbose
        fprintf(varargin{:});
    end
end
