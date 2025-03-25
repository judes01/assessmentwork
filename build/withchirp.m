f1 = 200; % Starting frequency of the chirp signal
f2 = 5000; % Ending frequency of the chirp signal
dur = 3; % Duration of the chirp signal in seconds
fs = 44100; % Sampling frequency in Hz
framesize = 2048; % Frame size for the spectrogram
overlap = 1024; % Overlap between frames for the spectrogram

% Generate the chirp signal using a custom function (ensure you have this function)
chirpsig = mychirp(f1, f2, dur, fs);

% Plot the waveform of the chirp signal
figure, wvtool(chirpsig)

% Plot the spectrograms
figure 

% Subplot 1: Spectrogram with rectangular window (default)
subplot(1, 2, 1)
spectrogram(chirpsig, framesize, overlap, framesize, fs)
title('Spectrogram with Rectangular Window')

% Subplot 2: Spectrogram with Hamming window
subplot(1, 2, 2)
spectrogram(chirpsig, hamming(framesize), overlap, framesize, fs)
title('Spectrogram with Hamming Window')

% Additional analysis: You might want to compare the frequency content of the chirp
% Optionally, you can plot the frequency response of the chirp signal:

figure
% Plot the frequency spectrum (magnitude) using FFT
N = length(chirpsig);
f = (0:N-1)*(fs/N); % Frequency axis
Y = abs(fft(chirpsig)); % Fourier Transform of the signal
plot(f(1:N/2), Y(1:N/2)) % Plot the magnitude spectrum for the first half (positive frequencies)
title('Frequency Spectrum of Chirp Signal')
xlabel('Frequency (Hz)')
ylabel('Magnitude')