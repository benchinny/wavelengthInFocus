function cmsg(str_, level_, verbosity_)
    % Just a small function for logging messages to the command window
    if nargin < 3
        verbosity_ = 1;
    end
    if nargin < 2
        level_ = 1;
    end
    if level_ == 5
        msg_ = 'CRITICAL';
    elseif level_ == 4
        msg_ = 'ERROR';
    elseif level_ == 3
        msg_ = 'WARNING';
    elseif level_ == 2
        msg_ = 'INFO';
    elseif level_ == 1
        msg_ = 'DEBUG';
    end
    if verbosity_ <= level_
        dt = datestr(now, 'yyyy-mm-dd HH:MM:SS.FFF');
        fprintf('%s : %s : %s\n', dt, msg_, str_)
    end
end