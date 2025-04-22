%% 1. Parameter Definition
% Room dimensions (meters)
room_length = 7.6;    % Length in x-direction
room_width  = 7.5;    % Width in y-direction
room_height = 4;      % Height in z-direction

% Sound source and receiver (listener) positions [x, y, z] in meters
src_pos = [3, 3, 1.5];   % Source position
rec_pos = [5, 5, 1.5];   % Receiver position

% Sampling parameters and speed of sound
fs = 44100;         % Sampling frequency in Hz
c  = 343;           % Speed of sound in m/s

% Maximum reflection order (i.e. how many reflections to calculate)
max_order = 50;

% Absorption coefficients for surfaces (values between 0 and 1)
% These represent the fraction of energy absorbed at each surface.
absorption_wall   = 0.1;   % Example for walls
absorption_floor  = 0.2;   % Example for floor
absorption_ceiling= 0.2;   % Example for ceiling
absorption_glass  = 0.15;  % Example for glass windows

% Impulse response duration (seconds)
ir_duration = 3.0;  
N = round(fs * ir_duration);    % Number of samples in the IR
h = zeros(N, 1);                % Preallocate the impulse response vector

%% 2. Image Source Method (ISM) and Impulse Response Construction
% Call the function to generate the impulse response
h = generate_IR(room_length, room_width, room_height, src_pos, rec_pos, fs, max_order, ...
                absorption_wall, absorption_floor, absorption_ceiling, absorption_glass, N, c);

%% 3. Plot the Impulse Response
time_axis = (0:N-1) / fs;  % Time vector in seconds
figure;
stem(time_axis, h, 'Marker', 'none');
xlabel('Time (s)');
ylabel('Amplitude');
title('Room Impulse Response using Image Source Method');
grid on;

%% 4. Using the Impulse Response for Convolution Reverb
% The generated impulse response (vector h) can now be used to apply convolution
% reverb to an audio signal.
% Example (uncomment and modify as needed):
%
%   [audio, fs_audio] = audioread('your_audio_file.wav');
%   audio_reverberated = conv(audio, h);
%   soundsc(audio_reverberated, fs);
function h = generate_IR(room_length, room_width, room_height, src_pos, rec_pos, fs, max_order, ...
                absorption_wall, absorption_floor, absorption_ceiling, absorption_glass, N, c)
    % Function to generate the impulse response using the Image Source Method (ISM)
    % Arguments:
    % - room dimensions: room_length, room_width, room_height
    % - source and receiver positions: src_pos, rec_pos
    % - sampling frequency: fs
    % - maximum reflection order: max_order
    % - absorption coefficients for various surfaces
    % - impulse response length: N
    % - speed of sound: c
    
    % Initialize impulse response vector
    h = zeros(N, 1);

    % Reflection coefficients for each surface
    reflection_coeff_wall = 1 - absorption_wall;
    reflection_coeff_floor = 1 - absorption_floor;
    reflection_coeff_ceiling = 1 - absorption_ceiling;
    reflection_coeff_glass = 1 - absorption_glass;

    % Loop over reflection orders in each spatial dimension (x, y, z)
    for nx = -max_order:max_order
        for ny = -max_order:max_order
            for nz = -max_order:max_order
                % Compute image source x-coordinate
                if mod(nx, 2) == 0
                    img_x = src_pos(1) + nx * room_length;
                else
                    img_x = (room_length - src_pos(1)) + nx * room_length;
                end

                % Compute image source y-coordinate
                if mod(ny, 2) == 0
                    img_y = src_pos(2) + ny * room_width;
                else
                    img_y = (room_width - src_pos(2)) + ny * room_width;
                end

                % Compute image source z-coordinate
                if mod(nz, 2) == 0
                    img_z = src_pos(3) + nz * room_height;
                else
                    img_z = (room_height - src_pos(3)) + nz * room_height;
                end

                % Compute the Euclidean distance from the image source to the receiver
                distance = sqrt((img_x - rec_pos(1))^2 + (img_y - rec_pos(2))^2 + (img_z - rec_pos(3))^2);

                % Compute the propagation delay in seconds and convert to a sample index
                time_delay = distance / c;
                sample_delay = round(time_delay * fs) + 1;  % +1 for MATLAB 1-indexing

                % Count the total number of reflections from all three dimensions
                num_reflections = abs(nx) + abs(ny) + abs(nz);

                % Compute the overall reflection coefficient for this image source
                % Reflection coefficient depends on the surface it reflects off of
                if mod(nx, 2) == 0
                    reflection_coeff = reflection_coeff_wall;  % Wall reflection
                elseif mod(ny, 2) == 0
                    reflection_coeff = reflection_coeff_wall;  % Wall reflection
                elseif mod(nz, 2) == 0
                    reflection_coeff = reflection_coeff_ceiling;  % Ceiling reflection
                else
                    reflection_coeff = reflection_coeff_floor;  % Floor reflection
                end

                % Attenuation due to spherical spreading (inverse of distance)
                attenuation = 1 / distance;

                % Add the contribution from this image source to the impulse response
                if sample_delay <= N
                    h(sample_delay) = h(sample_delay) + reflection_coeff^num_reflections * attenuation;
                end
            end
        end
    end
end
