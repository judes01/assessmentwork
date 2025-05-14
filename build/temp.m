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
% 5. Plots their frequency content using FFT

%% 1. SELECT AUDIO FILES

%clear; clc; close all% 


[filename1, pathname1] = uigetfile('*.wav', 'Select the FIRST audio file'); % 
if isequal(filename1, 0)
disp('User canceled first file selection.');
return;
end

[filename2, pathname2] = uigetfile('*.wav', 'Select the SECOND audio file');
if isequal(filename2, 0)
disp('User canceled second file selection.');
return;
end

% Read the audio files
[x1, fs1] = audioread(fullfile(pathname1, filename1));
[x2, fs2] = audioread(fullfile(pathname2, filename2));

% Check sample rate match
if fs1 ~= fs2
%error('Sample rates of the two audio files do not match.');
x2 = resample(x2, fs2, fs1);
end


%fs2 = fs1; % Use common sample rate

% Convert stereo to mono (if needed)
if size(x1, 2) == 2
x1 = mean(x1, 2);
end
if size(x2, 2) == 2
x2 = mean(x2, 2);
end

%% 2. GET FILTER SETTINGS FROM USER

filterType = input('Enter filter type (low or high): ', 's');
if ~ismember(filterType, {'low', 'high'})
error('Invalid filter type. Choose "low" or "high".');
end

fc = input('Enter cutoff frequency in Hz (e.g. 1000): ');
%fc = str2double(fc);
order = input('Enter filter order (e.g. 4): ');
ripple = input('Enter ripple in dB for Chebyshev filter (e.g. 1): ');
nyq = fs1 / 2;

Wn = fc / nyq; % Normalize cutoff frequency (0–1)

% Design Chebyshev Type I filter
[b, a] = cheby1(order, ripple, Wn, filterType);

%% 3. APPLY FILTER AND NORMALIZE

% Filter each signal
y1 = filter(b, a, x1);
y2 = filter(b, a, x2);

% Normalize (scale between -1 and 1)
y1 = y1 / max(abs(y1));
y2 = y2 / max(abs(y2));

% Save filtered files
audiowrite('filtered_file1.wav', y1, fs1);
audiowrite('filtered_file2.wav', y2, fs1);
disp('Both files filtered and saved.');

%% 4. PLOT MAGNITUDE SPECTRUM USING FFT

% Get signal lengths
N1 = length(y1);
N2 = length(y2);

% Apply FFT
Y1 = fft(y1);
Y2 = fft(y2);

% Get magnitude of positive frequencies only
Y1_mag = abs(Y1(1:N1/2+1));
Y2_mag = abs(Y2(1:N2/2+1));

% Frequency axis
f1 = (0:N1/2) * fs1 / N1;
f2 = (0:N2/2) * fs1 / N2;

% Normalize FFT magnitude
Y1_mag = Y1_mag / N1;
Y2_mag = Y2_mag / N2;

% Double middle bins to keep energy (except DC and Nyquist)
Y1_mag(2:end-1) = 2 * Y1_mag(2:end-1);
Y2_mag(2:end-1) = 2 * Y2_mag(2:end-1);

% Plot
figure;

subplot(2,1,1);
plot(f1, Y1_mag, 'b');
title(['Magnitude Spectrum: ', filename1]);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([0 nyq]);
grid on;

subplot(2,1,2);
plot(f2, Y2_mag, 'r');
title(['Magnitude Spectrum: ', filename2]);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([0 nyq]);
grid on;

drawnow; % Forces figure to render now

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
