function iFrequency = ifreqDeconv(signal,fs,smoothingWindow)
% calculates event weighted frequency of events that are separated by 0s
%
%signal = ciaDeconv(roi,:);
%signal = A;

if size(signal,1) > size(signal,2)
    signal = signal';
end

signal(signal>1) = 0; % eliminate artificial high values

iFrequency = nan(size(signal));

for i = 1:1:size(iFrequency,1)
    
    eventPos = find(signal(i,:) ~= 0);
    events = signal(i,eventPos);
    
    freques = fs./diff(eventPos);
    freques = [freques 0];                              % add 0 as the last datapoint because frequency can only be calculated between two events
    
    iFrequency(i,signal(i,:) ~= 0) = freques.*events;          % fill data in the correct positions, frequency is weighted by the event size
    iFrequency(i,1) = 0;                                  % add 0 as the first value for the fillmissing function
    iFrequency(i,:) = fillmissing(iFrequency(i,:),'previous');
    
    % divide data with the smalest "unitary" deconvolved value
    if max(signal(i,:)) > 0
        iFrequency(i,:) = iFrequency(i,:) / min(signal(i,(signal(i,:) > 0)));
       
    end
end

%iFrequency = smoothdata(iFrequency,2,'movmedian',3); % filtering high outlier events
iFrequency(iFrequency > 200) = 200;
iFrequency = smoothdata(iFrequency,2,'gaussian',smoothingWindow); % Note: the std for gaussian method is fixed 1/5th of the window 
end


%% Backup older working version (2020.06)
%{ 
function iFrequency = ifreqDeconv(signal,fs)
% calculates event weighted frequency of events that are separated by 0s
% 
%signal = ciaDeconv(roi,:);
%signal = A;

iFrequency = nan(size(signal));

if size(signal,1) > size(signal,2)
    signal = signal';
end

eventPos = find(signal ~= 0);
events = signal(eventPos);                          % first datapoint is excluded

freques = fs./diff(eventPos);
freques = [freques 0];                              % add 0 as the last datapoint because frequency can only be calculated between two events 

iFrequency(signal ~= 0) = freques.*events;          % fill data in the correct positions 
iFrequency(1) = 0;                                  % add 0 as the first value for the fillmissing function
iFrequency = fillmissing(iFrequency,'previous'); 

for i = 1:1:size(iFrequency,1)                      % divide data with the smapest "unitary" deconvolved value
    if max(signal(i,:)) > 0
        iFrequency(i,:) = iFrequency(i,:) / min(signal(i,(signal(i,:) > 0)));
    end
end

end

%}


