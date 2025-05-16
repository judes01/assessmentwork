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

%% ASSESSMENT SCRIPT — CHEBYSHEV FILTER AND FREQUENCY ANALYSIS

% This script:
% 1. Asks the user to select two audio files
% 2. Applies a Chebyshev Type I filter (low-pass or high-pass)
% 3. Normalizes the filtered audio
% 4. Saves the filtered versions
% 5. Plots both the magnitude and phase spectra


% Base functions used:
% audioread()   - Reads audio files
% audiowrite()  - Saves audio files
% filter()      - Applies a digital filter to audio
% normalize()   - Scales audio to fit within [-1, 1]
% fft()         - Converts audio from time domain to frequency domain

%% 1. SELECT AUDIO FILES

clear;    % Clear all variables from workspace
clc;      % Clear the command window

% Ask the user to choose the first audio file (.wav)
[filename1, pathname1] = uigetfile('*.wav', 'Select the FIRST audio file'); 

% If the user presses cancel, stop the script
if isequal(filename1, 0)
    disp('User canceled first file selection.');
    return;  % Exit the script
end

% Ask the user to choose the second audio file (.wav)
[filename2, pathname2] = uigetfile('*.wav', 'Select the SECOND audio file');

% If the user presses cancel, stop the script
if isequal(filename2, 0)
    disp('User canceled second file selection.');
    return;  % Exit the script
end

% Read the first audio file
[x1, fs1] = audioread(fullfile(pathname1, filename1));  % x1 is the audio data, fs1 is the sample rate

% Read the second audio file
[x2, fs2] = audioread(fullfile(pathname2, filename2));  % x2 is the audio data, fs2 is the sample rate

% Check if sample rates match, if not, resample x2 to match fs1
if fs1 ~= fs2
    x2 = resample(x2, fs2, fs1);  % Adjusts x2 so both files have same sample rate
end

% Convert stereo to mono by averaging left and right channels (if needed)
if size(x1, 2) == 2  % Check if x1 has 2 channels (stereo)
    x1 = mean(x1, 2);  % Convert to mono by averaging
end
if size(x2, 2) == 2  % Check if x2 has 2 channels (stereo)
    x2 = mean(x2, 2);  % Convert to mono by averaging
end

%% 2. GET FILTER SETTINGS FROM USER

% Ask the user whether they want a low-pass or high-pass filter
filterType = input('Enter filter type (low or high): ', 's');

% Validate input; if it's not "low" or "high", stop the script
if ~ismember(filterType, {'low', 'high'})
    error('Invalid filter type. Choose "low" or "high".');  % Show error and exit
end

% Ask the user for the cutoff frequency in Hz
fc = input('Enter cutoff frequency in Hz (e.g. 1000): ');

% Ask for filter order (higher = sharper filter)
order = input('Enter filter order (e.g. 4): ');

% Ask how much ripple (variation) is allowed in the passband, in dB
ripple = input('Enter ripple in dB for Chebyshev filter (e.g. 1): ');

% Calculate Nyquist frequency (half of sampling rate)
nyq = fs1 / 2;

% Convert cutoff frequency into a value between 0 and 1 (normalized)
Wn = fc / nyq;

% Design the Chebyshev Type I filter using user settings
[b, a] = cheby1(order, ripple, Wn, filterType);  % b = numerator, a = denominator coefficients

%% 3. APPLY FILTER AND NORMALIZE

% Apply the filter to the first audio signal
y1 = filter(b, a, x1);  % Filtered signal y1

% Apply the filter to the second audio signal
y2 = filter(b, a, x2);  % Filtered signal y2

% Normalize y1 so its loudest point is 1 or -1
y1 = y1 / max(abs(y1));  % Prevents clipping

% Normalize y2 so its loudest point is 1 or -1
y2 = y2 / max(abs(y2));  % Prevents clipping

% Save the filtered signals as new audio files
audiowrite('filtered_file1.wav', y1, fs1);  % Save y1

audiowrite('filtered_file2.wav', y2, fs1);  % Save y2

disp('Both files filtered and saved.');  % Let the user know it's done

%% 4. PLOT MAGNITUDE & PHASE SPECTRUM USING FFT

% Get the length of the first filtered signal
N1 = length(y1);

% Get the length of the second filtered signal
N2 = length(y2);

% Perform Fast Fourier Transform to move to frequency domain
Y1 = fft(y1);  % FFT of first signal
Y2 = fft(y2);  % FFT of second signal

% Get magnitude for only the positive frequencies
Y1_mag = abs(Y1(1:N1/2+1));
Y2_mag = abs(Y2(1:N2/2+1));

% Get phase angles for positive frequencies
Y1_phase = angle(Y1(1:N1/2+1));  % Phase of Y1
Y2_phase = angle(Y2(1:N2/2+1));  % Phase of Y2

% Create frequency axis for first signal
f1 = (0:N1/2) * fs1 / N1;

% Create frequency axis for second signal
f2 = (0:N2/2) * fs1 / N2;

% Normalize magnitude (to make comparison fair)
Y1_mag = Y1_mag / N1;
Y2_mag = Y2_mag / N2;

% Double all bins except the first and last to conserve signal energy
Y1_mag(2:end-1) = 2 * Y1_mag(2:end-1);
Y2_mag(2:end-1) = 2 * Y2_mag(2:end-1);

% Start a new figure window for plotting
figure;

% Plot magnitude spectrum for first file
subplot(2,2,1);
plot(f1, Y1_mag, 'b');
title(['Magnitude Spectrum: ', filename1]);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([0 nyq]);  % Limit x-axis to Nyquist frequency
grid on;

% Plot magnitude spectrum for second file
subplot(2,2,2);
plot(f2, Y2_mag, 'r');  
title(['Magnitude Spectrum: ', filename2]);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([0 nyq]);
grid on;

% Plot phase spectrum for first file
subplot(2,2,3);
plot(f1, unwrap(Y1_phase), 'b');  % Unwrap to avoid phase jumps
title(['Phase Spectrum: ', filename1]);
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
xlim([0 nyq]);
grid on;

% Plot phase spectrum for second file
subplot(2,2,4);
plot(f2, unwrap(Y2_phase), 'r');  % Unwrap phase
title(['Phase Spectrum: ', filename2]);
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
xlim([0 nyq]);
grid on;

drawnow;  % Render all plots immediately

%% EXPLANATION SECTION 

%% Background / application – why these techniques are necessary and where they can be
% found in the real world. 
% EXPLANATION (SHORT VERSION)
% - Chebyshev filters are sharp and allow ripple in the passband.
% - Normalization avoids distortion/clipping after filtering.
% - FFT shows which frequencies are present in each audio file


%% Theory – explanation of the theory and mathematics of how these processes achieve their
% results.



%% Features – What is your solution, what features are included




%% Code – discussion of your code implementation, how it works, how do the parameters effect
% the results.




%% Results – showing your solution solves the audio processing tasks



%% Limitations – what are the limitations, what bugs did you encounter, were you able to fix
% Error using newcode
% Sample rates of the two audio files do not match.
% 
% newcode
% Error using newcode
% Sample rates of the two audio files do not match.
% 
% newcode
% Error using newcode
% Sample rates of the two audio files do not match.
% 
% newcode
% Error using newcode
% Sample rates of the two audio files do not match.
% 
% newcode
% Enter filter type (low or high): low
% Enter cutoff frequency in Hz (e.g. 1000): 4000
% Enter filter order (e.g. 4): 100
% Enter ripple in dB for Chebyshev filter (e.g. 1): 1.8
% Unrecognized function or variable 'fs'.
% 
% Error in newcode (line 58)
% Wn = fc / (fs / 2); % Normalize cutoff frequency (0–1)
% 
% newcode
% Enter filter type (low or high): low
% Enter cutoff frequency in Hz (e.g. 1000): 4000
% Enter filter order (e.g. 4): 100
% Enter ripple in dB for Chebyshev filter (e.g. 1): 2
% Unrecognized function or variable 'fs'.
% 
% Error in newcode (line 59)
% Wn = fc / (fs / 2); % Normalize cutoff frequency (0–1)
% 
% newcode
% Enter filter type (low or high): low
% Enter cutoff frequency in Hz (e.g. 1000): 4000
% Enter filter order (e.g. 4): 100
% Enter ripple in dB for Chebyshev filter (e.g. 1): 2
% Error using lp2lp
% Expected input number 5, Wo, to be finite.
% 
% Error in lp2lp (line 61)
%     validateattributes(wo,{'numeric'},{'scalar','finite','real'},'lp2lp','Wo',5);
% 
% Error in cheby1 (line 108)
%    [a,b,c,d] = lp2lp(a,b,c,d,Wn);
% 
% Error in newcode (line 62)
% [b, a] = cheby1(order, ripple, Wn, filterType);
% 
% newcode
% Enter filter type (low or high): low
% Enter cutoff frequency in Hz (e.g. 1000): 4000
% Enter filter order (e.g. 4): 100
% Enter ripple in dB for Chebyshev filter (e.g. 1): 2
% newcode
% Enter filter type (low or high): low
% Enter cutoff frequency in Hz (e.g. 1000): 4000
% Enter filter order (e.g. 4): 100
% Enter ripple in dB for Chebyshev filter (e.g. 1): 2
% Warning: Matrix is close to singular or badly scaled. Results may be inaccurate. RCOND =  1.474178e-19. 
% > In bilinear (line 168)
% In cheby1 (line 122)
% In newcode (line 63)
% 
% Warning: Matrix is close to singular or badly scaled. Results may be inaccurate. RCOND =  1.474178e-19. 
% > In bilinear (line 169)
% In cheby1 (line 122)
% In newcode (line 63)
% 
% Warning: Matrix is close to singular or badly scaled. Results may be inaccurate. RCOND =  1.474178e-19. 
% > In bilinear (line 170)
% In cheby1 (line 122)
% In newcode (line 63)
% 
% Warning: Matrix is close to singular or badly scaled. Results may be inaccurate. RCOND =  1.474178e-19. 
% > In bilinear (line 171)
% In cheby1 (line 122)
% In newcode (line 63)
% 
% Unrecognized function or variable 'fs'.
% 
% Error in newcode (line 76)
% audiowrite('filtered_file1.wav', y1, fs);
% 
% newcode
% Enter filter type (low or high): low
% Enter cutoff frequency in Hz (e.g. 1000): 4000
% Enter filter order (e.g. 4): 100
% Enter ripple in dB for Chebyshev filter (e.g. 1): 2
% Unrecognized function or variable 'fs'.
% 
% Error in newcode (line 77)
% audiowrite('filtered_file1.wav', y1, fs);
% 
% newcode
% Enter filter type (low or high): low
% Enter cutoff frequency in Hz (e.g. 1000): 4000
% Enter filter order (e.g. 4): 100
% Enter ripple in dB for Chebyshev filter (e.g. 1): 2
% Warning: Matrix is close to singular or badly scaled. Results may be inaccurate. RCOND =  1.474178e-19. 
% > In bilinear (line 168)
% In cheby1 (line 122)
% In newcode (line 64)
% 
% Warning: Matrix is close to singular or badly scaled. Results may be inaccurate. RCOND =  1.474178e-19. 
% > In bilinear (line 169)
% In cheby1 (line 122)
% In newcode (line 64)
% 
% Warning: Matrix is close to singular or badly scaled. Results may be inaccurate. RCOND =  1.474178e-19. 
% > In bilinear (line 170)
% In cheby1 (line 122)
% In newcode (line 64)
% 
% Warning: Matrix is close to singular or badly scaled. Results may be inaccurate. RCOND =  1.474178e-19. 
% > In bilinear (line 171)
% In cheby1 (line 122)
% In newcode (line 64)
% 
% Both files filtered and saved.
% Warning: Integer operands are required for colon operator when used as index. 
% > In newcode (line 92)
% 
% newcode
% Enter filter type (low or high): low
% Enter cutoff frequency in Hz (e.g. 1000): 4000
% Enter filter order (e.g. 4): 6
% Enter ripple in dB for Chebyshev filter (e.g. 1): 2
% Both files filtered and saved.
% Warning: Integer operands are required for colon operator when used as index. 
% > In newcode (line 92)
