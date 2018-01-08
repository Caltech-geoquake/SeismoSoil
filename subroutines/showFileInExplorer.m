function showFileInExplorer(full_file_name)
%
% Opens Explorer window (win) or Finder window (mac) that points to the
% location speficied by full_file_name. For Windows, the specific file
% within the Explorer window can be highlighted.
%
% (c) 2018/1/7, Jian Shi

if ispc()
    command_Windows = sprintf('explorer /select,%s',full_file_name);
    system(command_Windows);
else
    cdir = fileparts(full_file_name);
    cdir = strrep(cdir,'(','\('); % replace ( with \(
    cdir = strrep(cdir,')','\)'); % replace ) with \)
    cdir = strrep(cdir,' ','\ '); % add a \ before the spaces
    command_Mac = sprintf('open %s',cdir);
    system(command_Mac);
end

end

