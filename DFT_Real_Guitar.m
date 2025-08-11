%DFT_Real_Guitar - This MATLAB script processes the audio signals of the
%                  open E, A, and D guitar strings. It loads the audio
%                  data, trims it to a 1-second segment, and applies the
%                  Hamming window. The DFT of the signals is computed with
%                  and without windowing, and visualized via plots. The
%                  plots help examine the frequency content of each
%                  string's sound
%                  


%% Loading in Audio files

Fs = 4410; % sets the sampling frequency
T = 1 / Fs; % sets the sampling interval

% Loading audio files into MATLAB
x_E = audioread('E_String wav file.wav');
x_A = audioread('A_String Wav File.wav');
x_D = audioread('D String Wav File.wav');

% converting audio vectors into row vectors
x_E = x_E.';
x_A = x_A.';
x_D = x_D.';

% obtaining length of audio vectors
N_D = length(x_D);
N_A = length(x_A);
N_E = length(x_E);

n_D = 0:(N_D -1);
n_A = 0:(N_A -1);
n_E = 0:(N_E -1);

t_D = n_D * T;
t_A = n_A * T;
t_E = n_E * T;

figure(1)
subplot(3,1,1)
plot(t_E, x_E, 'k');
xlabel('Time (s)'); ylabel('Signal x_E[n]');
title('Open Low E string');

subplot(3,1,2)
plot(t_A, x_A, 'r');
xlabel('Time (s)'); ylabel('Signal x_A[n]');
title('Open A string');

subplot(3,1,3)
plot(t_D, x_D);
xlabel('Time (s)'); ylabel('Signal x_D[n]')
title('Open D string')

% saving the 5th second of the notes (between t = 4 and t = 5)
x_E = x_E((t_E >= 4) & (t_E <= 5));
x_D = x_D((t_D >= 4) & (t_D <= 5));
x_A = x_A((t_A >= 4) & (t_A <= 5));

% we need to recalculate the length of the truncated signal,
% however we can calculate the new length after truncation to be 4411
% samples

N = 4411;

% plotting the truncated signals 
figure(2)
subplot(3,1,1)
plot(t_E((t_E >= 4) & (t_E <= 5)), x_E, 'k');
xlabel('Time (s)'); ylabel('Signal x_E[n]');
title('Open Low E string (Truncated)');

subplot(3,1,2)
plot(t_A((t_D >= 4) & (t_D <= 5)), x_A, 'r');
xlabel('Time (s)'); ylabel('Signal x_A[n]');
title('Open A string (Truncated)');

subplot(3,1,3)
plot(t_D((t_A >= 4) & (t_A <= 5)), x_D);
xlabel('Time (s)'); ylabel('Signal x_D[n]')
title('Open D string (Truncated)')

%% Computing the DFT of the audio signals

w = hamming(N);  % hamming window for E string

% obtain hamming vector as a row vector
w = w.';

% creates windowed singal 
xw_E = x_E .* w;
xw_A = x_A .* w;
xw_D = x_D .* w;

Nzp = 32768; % amount of zero-padding to be used when computing DFT

% Computing the DFT of the signals with and without windowing
X_A = fftshift(fft(x_A, Nzp));
XW_A = fftshift(fft(xw_A, Nzp));

X_E = fftshift(fft(x_E, Nzp));
XW_E = fftshift(fft(xw_E, Nzp));

X_D = fftshift(fft(x_D, Nzp));
XW_D = fftshift(fft(xw_D, Nzp));

freq = ((-Nzp/2):((Nzp/2)-1)) * (Fs/Nzp);  % creates frequency points for plot in Hz

%% Plotting the DFTs 
figure(3)
subplot(2,1,1)
plot(freq, 20*log10(abs(X_E)), 'k');
grid on 
xlabel('Frequency (Hz)');
ylabel('|X(\Omega)| (dB)');
title('DFT of Open Low E string without windowing');
axis([0 1000 -80 50])

subplot(2,1,2)
plot(freq, 20*log10(abs(XW_E)), 'k');
grid on 
xlabel('Frequency (Hz)');
ylabel('|X(\Omega)| (dB)');
title('DFT of Open Low E string with windowing');
axis([0 1000 -80 50])

figure(4)
subplot(2,1,1)
plot(freq, 20*log10(abs(X_A)), 'r');
grid on 
xlabel('Frequency (Hz)'); ylabel('|X(\Omega)| (dB)');
title('DFT of Open A string without windowing');
axis([0 1000 -80 50])

subplot(2,1,2)
plot(freq, 20*log10(abs(XW_A)), 'r');
grid on 
xlabel('Frequency (Hz)'); ylabel('|X(\Omega)| (dB)');
title('DFT of Open A string with windowing');
axis([0 1000 -80 50])

figure(5)
subplot(2,1,1)
plot(freq, 20*log10(abs(X_D)));
grid on 
xlabel('Frequency (Hz)'); ylabel('|X(\Omega)| (dB)');
title('DFT of Open D string without windowing');
axis([0 1000 -80 50])

subplot(2,1,2)
plot(freq, 20*log10(abs(XW_D)));
grid on 
xlabel('Frequency (Hz)'); ylabel('|X(\Omega)| (dB)');
title('DFT of Open D string with windowing');
axis([0 1000 -80 50])