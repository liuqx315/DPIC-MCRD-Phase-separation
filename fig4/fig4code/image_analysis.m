
clear all
files = dir( "*.tif");
filenames = {files.name};
Tmax = numel(filenames);
for i = 1:17

    spectral2D(filenames{i},'b=1,p=1');
    
end
close all