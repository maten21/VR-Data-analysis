% Simple way to extract indexes of certain trial types from the trial type array
% Use "t" as index



trials = find(sData.trials.trialTypeArrayStim == sData.trials.trialTypesMeta(t).trialTypeIndicator);
nTrials = numel(trials);
