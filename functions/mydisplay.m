function mydisplay(verbose, varargin)
    if nargin == 3
        text = strcat(varargin{1}, " = ");
        table = varargin{2};
    elseif nargin == 2
        text = "";
        table = varargin{1};
    else
        error("Incompatible number of optional arguments.")
    end

    if verbose
        fprintf(text);
        display(table);
    end
end