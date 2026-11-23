function out = custom_method_template(c,cfg,params) %#ok<INUSD>
% Copy/rename this function and replace the intentional error with your method.
% Register: cfg.methods(end+1)=struct('name','my_method', ...
%    'run',@custom_method_template,'params',struct());
% c contains ONE muscle: .envelope, .y, .energy, .valid, .dynamicValid,
% .t, .fs, .mvc; .raw/.rawTime/.rawFs/.rawValid allow another demodulator.
% Preserve the output time grid and all invalid intervals. Map amplitudes to
% the common reference-envelope domain before comparing physical powers.
out=tbench.emptyOutput(c); %#ok<NASGU>
error('tbench:Template','Implement this template before registering it.');
% For static output: frequency, powerMean, detected; leave dynamic fields empty.
% For dynamic output: capabilities.dynamic=true, capabilities.waveform=true;
% signal/fundamental (voltage), active (logical), frequencyTrack (Hz), valid.
% onProbability is optional and must not be called a calibrated event probability.
end
