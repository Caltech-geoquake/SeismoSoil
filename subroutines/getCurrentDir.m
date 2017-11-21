function currentDir = getCurrentDir
%if isdeployed % Stand-alone mode.
    if strcmpi(computer,'pcwin64')
        [status,result] = system('path');
    elseif strcmpi(computer,'pcwin')
        [status,result] = system('path');
    elseif strcmpi(computer,'maci64')
        [status,result] = system('echo $PATH');
    end
    if status == 0
        currentDir = char(regexpi(result, 'Path=(.*?);', 'tokens', 'once'));
    end
%else % MATLAB mode.
%    currentDir = pwd;
end