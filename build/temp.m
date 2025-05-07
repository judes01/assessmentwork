%%assessment work

audioread() % Reads an audio file and returns the audio data and sample rate.
audiowrite() % Writes audio data to a file with the specified format and sample rate.
conv() % Performs convolution of two signals or arrays. It computes the sum of the product of two sequences over a range of time.
filter() % Applies a filter to a signal (e.g., low-pass, high-pass) by using the given filter coefficients.
normalize() % Adjusts the amplitude of an audio signal to ensure it stays within a specified range (typically [-1, 1]).

 % i have to use the functions above for my code 
%applyng a filter
y = filter(b, a, x);
fs = 44100;
fc = 100;
[b, a] = butter(2, fc/(fs2));

%task
https://www.dspguide.com/pdfbook.htm
% simple_filter_cheby.m
% This script applies a Chebyshev Type I filter (low-pass or high-pass) to an audio file.
% It reads the file, filters the signal, normalizes the output, and saves the result.

% Prompt the user to select an audio file from their computer
[filename, pathname] = uigetfile('*.wav', 'Select an audio file to filter');

% Check if the user canceled the file selection
if isequal(filename,0)
    disp('User canceled the file selection.');
    return;  % Stop the script
end

% Read the selected audio file. 'x' is the audio data, 'fs' is the sample rate
[x, fs] = audioread(fullfile(pathname, filename)); % audioread reads the audio data

% If the audio is stereo (2 channels), separate it into left and right channels
if size(x,2) == 2
    left = x(:,1);   % Left channel
    right = x(:,2);  % Right channel
else
    left = x;        % If it's mono, use it as the left channel
    right = [];      % No right channel
end

% Ask the user to enter the type of filter they want to apply
filterType = input('Enter filter type (low or high): ', 's');

% Check that the user entered either "low" or "high"
if ismember(filterType, {'low', 'high'})
    % If the input is valid, continue
else
    % If the input is invalid, stop the script and show an error
    error('Invalid filter type. Choose "low" or "high".');
end

% Ask the user for the cutoff frequency of the filter (in Hz)
fc = input('Enter cutoff frequency in Hz (e.g. 1000): ');

% Ask the user for the filter order (higher means steeper but more complex)
order = input('Enter filter order (e.g. 4): ');

% Ask the user how much ripple is allowed in the passband (in dB)
% This is a unique feature of Chebyshev Type I filters
ripple = input('Enter ripple in dB for Chebyshev filter (e.g. 1): ');

% Convert the cutoff frequency to a value between 0 and 1 (required by MATLAB)
% This is done by dividing by half the sampling rate (Nyquist frequency)
Wn = fc / (fs / 2);

% Design the Chebyshev Type I filter using the user inputs
% 'b' and 'a' are the filter coefficients used to apply the filter
[b, a] = cheby1(order, ripple, Wn, filterType);

% Apply the filter to the audio signal using the filter coefficients
% This reduces or removes unwanted frequencies from the audio
y = filter(b, a, x);  % 'filter' applies the designed filter to the signal

% Normalize the filtered signal so that its values stay within the range [-1, 1]
% This prevents distortion and ensures proper playback volume
y = y / max(abs(y));  % Scale signal to maximum amplitude of 1

% Save the filtered audio signal to a new file called "filteredsig.wav"
% The file will be saved in the current MATLAB working directory
audiowrite('filteredsig.wav', y, fs);  % Write the output using the same sample rate

% Tell the user the process is finished and the file was saved
disp('Filtered audio saved as "filteredsig.wav".');

%% explainations 
%Chebyshev Type I filters are useful when you want a sharper cutoff than Butterworth filters, but can tolerate small ripples in the passband.
%Normalization ensures safe audio levels after processing.
%User interaction makes the script easy for beginners — no need to hard-code anything.
