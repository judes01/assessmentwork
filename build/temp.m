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
[b, a] = cheby1(2, fc/(fs2));

%task
https://www.dspguide.com/pdfbook.htm
%% simple_filter_cheby.m
% This script applies a Chebyshev Type I filter (low-pass or high-pass) to an audio file.
% It reads the file, filters the signal, normalizes the output, and saves the result.


%%1. Select Two Audio Files for Comparison
[filename1, pathname1] = uigetfile('*.wav', 'Select the FIRST audio file');  % Opens a dialog to select the first audio file
if isequal(filename1, 0)  % Checks if the user canceled the file selection
    disp('User canceled first file selection.');  % Displays a message if the user canceled
    return;  % Exits the script
end

[filename2, pathname2] = uigetfile('*.wav', 'Select the SECOND audio file');  
if isequal(filename2, 0)  
    disp('User canceled second file selection.');  
    return;  
end

% Read the audio files
[x1, fs1] = audioread(fullfile(pathname1, filename1));  % Reads the first audio file and its sample rate
[x2, fs2] = audioread(fullfile(pathname2, filename2)); 

% Ensure both files have the same sample rate
if fs1 ~= fs2  % Checks if the sample rates of the two files are different
    error('Sample rates of the two audio files do not match.');  % Throws an error if the sample rates don't match
end
fs = fs1;  % Uses the common sample rate for further processing

% Convert stereo to mono if necessary
if size(x1, 2) == 2  % Checks if the first audio file is stereo
    x1 = mean(x1, 2);  % Converts stereo to mono by averaging the channels
end
if size(x2, 2) == 2  
    x2 = mean(x2, 2); 
end

%%2. Ask user for filter settings
filterType = input('Enter filter type (low or high): ', 's');  % Prompts user for filter type (low or high)
if ~ismember(filterType, {'low', 'high'})  % Checks if the entered filter type is valid
    error('Invalid filter type. Choose "low" or "high".');  % Throws an error if the filter type is invalid
end

fc = input('Enter cutoff frequency in Hz (e.g. 1000): ');  % Prompts user for cutoff frequency
order = input('Enter filter order (e.g. 4): ');  % Prompts user for filter order
ripple = input('Enter ripple in dB for Chebyshev filter (e.g. 1): ');  % Prompts user for ripple value
Wn = fc / (fs / 2);  % Normalizes the cutoff frequency

% Design the filter
[b, a] = cheby1(order, ripple, Wn, filterType);  % Designs a Chebyshev filter based on user input

% Filter and normalize
y1 = filter(b, a, x1);  % Filters the first audio signal
y2 = filter(b, a, x2);  % Filters the second audio signal
y1 = y1 / max(abs(y1));  % Normalizes the first filtered signal
y2 = y2 / max(abs(y2));  % Normalizes the second filtered signal

% Save filtered files (optional, but useful for validation)
audiowrite('filtered_file1.wav', y1, fs);  % Saves the first filtered file to disk
audiowrite('filtered_file2.wav', y2, fs);  % Saves the second filtered file to disk
disp('Both files filtered and saved.');  % Displays a confirmation message

%%3. FFT and Magnitude Spectrum

% Get the length of each audio signal (y1 and y2)
N1 = length(y1);  % Gets the length of the first filtered signal
N2 = length(y2);  % Gets the length of the second filtered signal

% Perform the Fast Fourier Transform (FFT) on both audio signals
Y1 = fft(y1);  % Computes the dft of the first filtered signal
Y2 = fft(y2);  % Computes the dft of the second filtered signal

% Compute the magnitude of the FFT for both signals (get rid of negative frequencies)
Y1_mag = abs(Y1(1:N1/2+1));  % Get the magnitude of the positive frequencies of the first signal
Y2_mag = abs(Y2(1:N2/2+1));  % Get the magnitude of the positive frequencies of the second signal

% Generate the frequency axis (the range of frequencies we're interested in)
f1 = (0:N1/2) * fs / N1;  % Create the frequency axis for the first signal
f2 = (0:N2/2) * fs / N2;  % Create the frequency axis for the second signal

% Normalize the magnitude spectra (so the values are between 0 and 1)
Y1_mag = Y1_mag / N1;  % Normalize the magnitude of the first signal
Y2_mag = Y2_mag / N2;  % Normalize the magnitude of the second signal

% Scale the middle part of the FFT to account for the energy in the negative frequencies
Y1_mag(2:end-1) = 2 * Y1_mag(2:end-1);  % Double the middle part of the magnitude of the first signal (exclude the first and last values)
Y2_mag(2:end-1) = 2 * Y2_mag(2:end-1);  % Double the middle part of the magnitude of the second signal (exclude the first and last values)

% Create the figure to plot the graphs
figure;  % Start a new figure window for the plots

% Plot the magnitude spectrum of the first signal (filtered)
subplot(2,1,1);  % Create the first subplot (top plot in a 2-row layout)
plot(f1, Y1_mag, 'b');  % Plot the frequency vs magnitude of the first signal in blue
title(['Magnitude Spectrum (Filtered): ', filename1]);  % Set the title with the filename of the first signal
xlabel('Frequency (Hz)');  % Label the x-axis as frequency (in Hz)
ylabel('Magnitude');  % Label the y-axis as magnitude
xlim([0 fs/2]);  % Limit the x-axis to the range from 0 to half the sample rate (Nyquist frequency)

% Plot the magnitude spectrum of the second signal (filtered)
subplot(2,1,2);  % Create the second subplot (bottom plot in a 2-row layout)
plot(f2, Y2_mag, 'r');  % Plot the frequency vs magnitude of the second signal in red
title(['Magnitude Spectrum (Filtered): ', filename2]);  % Set the title with the filename of the second signal
xlabel('Frequency (Hz)');  % Label the x-axis as frequency (in Hz)
ylabel('Magnitude');  % Label the y-axis as magnitude
xlim([0 fs/2]);  % Limit the x-axis to the range from 0 to half the sample rate (Nyquist frequency)

%% explainations 
%Chebyshev Type I filters are useful when you want a sharper cutoff than Butterworth filters, but can tolerate small ripples in the passband.
%Normalization ensures safe audio levels after processing.
%User interaction makes the script easy for beginners — no need to hard-code anything.
