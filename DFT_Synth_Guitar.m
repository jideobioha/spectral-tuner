function [x, X] = DFT_Synth_Guitar(fo)
%DFT_Synth_Guitar - The purpose of this function is generate a synthetic
%                   audio signal similar to a guitar string given some base
%                   parameters for the signal. The signal is then windowed
%                   and its DFT is plotted, along with another plot of the 
%                   signal's DFT without windowing. 
% 
%
%   USAGE: [x,X] = DFT_Synth_Guitar(fo)
%   
%   inputs: fo - The fundamental frequency of the sinusoid in Hz
%
%   outputs:  This function creates plots of the DFT of the synthetic
%             signal with and without windowing. This function also returns
%             the synthetic signal x and the vector containing its DFT, X. 

Fs = 4410; %sets the sampling frequency
T = 1/Fs; %sets the samplimg period
n = 0:4410; % samples collected for a 1 second period with a frequency of 4410Hz
t = n * T;

phase = 2*pi*rand(1, 10); % creates a vector of random phase shifts

% create rw vector of amplitudes
Amp = [0.07 0.025 0.013 0.007 0.003 0.002 0.0004 0.0003 0.0002 0.0001];

x = 0; % initializing the variable that will store the signal

for k = 1:length(phase)
    current_term = Amp(k) * cos((k*2*pi*fo*t) + phase(k));
    x = x + current_term;
end

Nzp = 32768; % amount of zero-padding to be used when computing DFT

% Obtaining hamming window vector
w = hamming(length(x));
w = w.';

%creating windowed signal 
xw = x .* w;

%computing DFT of signals
X = fftshift(fft(x, Nzp));
XW = fftshift(fft(xw, Nzp));

%% Plotting 

freq = ((-Nzp/2):((Nzp/2)-1)) * (Fs/Nzp);  % creates frequency points for plot in Hz

figure;
plot(t, x); title('Synthetic Signal'); xlabel('Time (seconds)'); ylabel('Amplitude')

figure;
subplot(2,1,1)
plot(freq, 20*log10(abs(X)));
grid on 
xlabel('Frequency (Hz)'); ylabel('|X(\Omega)| (dB)');
title('DFT of Synthetic Guitar string without windowing');
axis([0 1000 -80 50])

subplot(2,1,2)
plot(freq, 20*log10(abs(XW)));
grid on 
xlabel('Frequency (Hz)'); ylabel('|X(\Omega)| (dB)');
title('DFT of Synthetic Guitar string with windowing');
axis([0 1000 -80 50])