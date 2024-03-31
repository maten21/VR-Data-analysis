function nFrames = numTiffDirectories(tiffRef)


    if ischar(tiffRef) && isfile(tiffRef)
        [~, ~, ext] = fileparts(tiffRef);
        if strcmp(ext, '.tif') || strcmp(ext, '.tiff')
            tiffObj = Tiff(tiffRef);
        else
            error('First input must be the path to an existing tiff file.')
        end
    elseif isa(tiffRef, 'Tiff')
        tiffObj = tiffRef;
    else
        error('First input must be a Tiff object or the path to a tiff file')
    end



try
    tiffObj.setDirectory(2);  
    offset2 = tiffObj.getTag(273); % get offset tag

    tiffObj.setDirectory(1);
    offset1 = tiffObj.getTag(273);

    frameSize = offset2 - offset1;

    rawImageSize = tiffObj.getTag(279); % Strip Byte Counts Field % tiffObj.getTag(256) * tiffObj.getTag(257) * tiffObj.getTag(258)/8; % imageResXY, bitPerPixel
    frameHeaderSize = frameSize - rawImageSize;
    mainHeader = offset1 - frameHeaderSize;

    d = dir(tiffObj.FileName);

    nFrames = double(d.bytes - mainHeader) / double(frameSize);
    
    if nFrames ~= round(nFrames) % in some files it looks like the start of the frameHeader is indicated by the offset
        nFrames = double(d.bytes - offset1) / double(frameSize);
    end

catch
    % if tiffObj.setDirectory() can not be set to 2 
    nFrames = 1;

end

if nFrames ~= round(nFrames)
    msgbox('The resulted nFrames is not a round number, the report error.')
    
    nFrames = floor(nFrames);

end

end



