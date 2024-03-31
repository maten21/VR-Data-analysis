function [hFig] = tif2png

[file, path] = uigetfile('*.tif','','MultiSelect','off' );

filePath = fullfile(path,file);

tiffObj = Tiff(filePath);
im = read(tiffObj);


hFig = figure('Position',[10, -600, size(im,2), size(im,2)],'Visible','off');
imshow(im)

[path, file, ext] = fileparts(filePath);
exportgraphics(hFig, fullfile(path, [file '3.png'])); 
end
