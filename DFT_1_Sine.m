function DFT_1_Sine(f1, N, Nzp)
%DFT_1_Sine - The purpose of this function is to compute the
%             DFT of a sinusoid with a frequency f1 in Hz with 
%             N samples taken from the sinusoid and Nzp samples
%             with zero padding. This function produces 2 plots. 
%             one plot representing DFT computed without using 
%             a window, and one plot representing the DFT computed 
%             using the Hamming Window
%
%   USAGE: DFT_1_Sine(f1, N, Nzp)
%   
% inputs: f1 - The frequency of the sinusoid in Hz
%           N  - The number the samples taken
%          Nzp - The total number of samples with zero padding
%  
%   outputs: This function outputs 2 plots of the DFT of a sinusoid
%             one plot with the DFT created without the Hamming Window
%             and one plot created with the Hamming Window

Fs = 4410; % Fs represents the Sampling Frequency
DT_f =((2 * pi) / (Fs)) * f1; % computes the DT frequency
n = 0:(N-1); % creates N DT sample index points
x = 1 * cos(DT_f * n); % Creates a row vector having N samples of a sinusoid with frequency of f1 Hz
X = fftshift(fft(x, Nzp)); % computes the DFT of DT signal x with zero-padding up to Nzp samples.
w = hamming(N); %creates Hamming Window vector as a column vector
w = w'; % transposes Hamming Window vector so as to obtain a row vector
xw = x .* w; % uses Hamming Window to create windowed signal 
XW = fftshift(fft(xw, Nzp));  % computes the DFT of windowed DT signal xw with zero-padding

T = 1/Fs; % set sampling period
t = n * T;
% plot signal
figure; plot(t,x); xlabel('t (seconds)'); ylabel('x(t)'); title('Input sinusoid');

freq = ((-Nzp/2):((Nzp/2)-1)) * (Fs/Nzp);  % creates frequency points for plot in Hz
nonZero_Frequencies = find(freq > 0);
%% Plots 

%Plotting the DFT of x without the Hamming Window
figure 
subplot(2,1,1)
plot(freq(nonZero_Frequencies), 20*log10(abs(X(nonZero_Frequencies))));
grid on 
xlabel('Frequency (Hz)')
ylabel('|X(\Omega)| (dB)')
title('DFT of Sinusoid computed without Hamming Window')

% Plotting the DFT of xw (with Hamming Window used)
subplot(2,1,2)
plot(freq(nonZero_Frequencies), 20*log10(abs(XW(nonZero_Frequencies))));
grid on 
xlabel('Frequency (Hz)')
ylabel('|X(\Omega)| (dB)')
title('DFT of Windowed Sinusoid Input Signal')

% Plotting the Hamming window
figure
plot(n, w)
grid on
xlabel('Time (s)')
ylabel('Hamming Window w')
title('Hamming Window Vector')

end