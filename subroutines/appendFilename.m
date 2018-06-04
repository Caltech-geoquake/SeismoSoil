function new_filename = appendFilename(original_filename, tail)
% Append a "tail" to the file name (before the file extension, after the
% original file name (without extension).

    [dir, name, ext] = fileparts(original_filename);
    new_filename = fullfile(dir, strcat(name, tail, ext));

end
