function openFolder(full_path,filename)
%
% Opens Explorer window (win) or Finder window (mac) that points to the
% location speficied by full_path. 
% 
% For Windows, the specified file (specified by the second input argument)
% within the Explorer window can be highlighted. For macOS, the second
% input argument is always ignored, because Terminal and Finder does not
% support such file highlighting.
%
% (c) 2018/1/7, Jian Shi

if ispc()
    if nargin == 2
        command_Windows = sprintf('explorer /select,%s',fullfile(full_path,filename));
    else
        command_Windows = sprintf('explorer %s',full_path);
    end
    system(command_Windows);
else
    cdir = fileparts(full_path);
    cdir = strrep(cdir,'(','\('); % replace ( with \(
    cdir = strrep(cdir,')','\)'); % replace ) with \)
    cdir = strrep(cdir,' ','\ '); % add a \ before the spaces
    command_Mac = sprintf('open %s',cdir);
    system(command_Mac);
end

end

