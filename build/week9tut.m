% Parameters for the chirp signal
Fs = 10000;        % Sampling frequency
A = 0.8;           % Amplitude of the chirp signal
Ts = 1/Fs;         % Sampling period
dur = 1.5;         % Duration of the chirp signal (in seconds)
f1 = 100;          % Starting frequency of the chirp (in Hz)
f2 = 5000;         % Ending frequency of the chirp (in Hz)

% Time vector
t = 0:Ts:dur;

% Instantaneous phase (theta) for a quadratic FM chirp
Theta = 2*pi*(f1 + (f2 - f1)/(2*dur) * t.^2); 

% Chirp signal using sin function
chirpsig = A * sin(Theta);

% Plot the waveform of the chirp signal
figure;
plot(t, chirpsig);
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Chirp Signal Waveform');

% Save the chirp signal to a WAV file
audiowrite('chirp_signal.wav', chirpsig, Fs);
disp('Chirp signal saved to chirp_signal.wav');

% Function to generate a linear FM chirp signal
function signal = mychirp(f1, f2, dur, fs)
    % MYCHIRP Generate a linear-FM chirp signal
    %
    % Usage: xx = mychirp(f1, f2, dur, fs)
    %
    % f1 = starting frequency (Hz)
    % f2 = ending frequency (Hz)
    % dur = total duration of the chirp signal (seconds)
    % fs = sampling frequency (Hz)

    % Check if fs is provided; if not, set default fs = 8000
    if nargin < 4
        fs = 8000; % default sampling frequency
    end
    
    % Time vector
    ts = 1 / fs; % Sampling period
    t = 0:ts:dur-ts; % Time vector
    
    % Calculate the frequency slope
    a = (f2 - f1) / (2 * dur); % Slope of the chirp
    
    % Generate the instantaneous phase (theta)
    theta = 2 * pi * (f1 * t + a * t.^2);
    
    % Generate the chirp signal by taking the real part of the complex exponential
    signal = real(exp(1i * theta));
end

