function iFrequency = ifreq(signal,fs)
% returns the instantaneous frequency of a discrete signal (containing 0 and 1) with fs sampling rate
% 
%signal = sData.behavior.signals.lickEvents;
%

iFrequency = nan(size(signal));

if size(signal,1) > size(signal,2)
    signal = signal';
end

events = find(signal == 1);
freques = fs./diff(events);
freques = [freques 0];

iFrequency(signal == 1) = freques;
iFrequency(1) = 0;
iFrequency = fillmissing(iFrequency,'previous'); 


end

