%% A personal function for saving figures into specified formats (such as colour eps).

% f is the handle of the figure to save
% filename is the name to give the file without the file extension
% extension is the file format in which to save the figure
function savefigas(f,filename,extension)

% Check if the chosen extension requires a particular driver
switch extension
    case 'eps'
        driver = 'epsc';
    otherwise
        driver = '';
end

% Build the full filename, including extension
filename = strcat(filename,'.',extension);

% Call the saveas function with the appropriate number of arguments
if isempty(driver)
    saveas(f,filename);
else
    saveas(f,filename,driver);
end

end