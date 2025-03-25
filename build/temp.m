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
