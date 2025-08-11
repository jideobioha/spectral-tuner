function DFT_2_Sine(f1, f2, a, N, Nzp)
%DFT_2_Sine - The purpose of this function is to compute the
%             DFT of 2 sinusoids with frequencies f1 and f2 in Hz with one
%             of the sinusoids having an amplitude of a with N sample
%             taken from the sinusoids and Nzp samples with zero padding
%             This function produces 2 plots. One plot representing DFT
%             computed without using a window, and one plot representing
%             representing the DFT computed using the Hamming window
%
%   USAGE: DFT_2_Sine(f1, f2, a, N, Nzp)
%   
% inputs: f1 - The frequency of one sinusoid in Hz
%         f2 - The frequency of the other sinusoid in Hz
%          a - Amplitude of the other sinusoid in volts
%           N  - The number the samples taken
%          Nzp - The total number of samples with zero padding
%  
%   outputs: This function outputs 2 plots of the DFT of a sinusoid
%             one plot with the DFT created without the Hamming Window
%             and one plot created with the Hamming Window

Fs = 4410; % Fs represents the Sampling Frequency
DT_f1 = ((2 * pi) / (Fs)) * f1; % computes the DT frequency for the 1st sinusoid
DT_f2 = ((2 * pi) / (Fs)) * f2; % comptues the DT frequency for the 2nd sinusoid
n = 0:(N-1); % creates N DT sample index points
x = (1 * cos(DT_f1 * n)) + (a * cos(DT_f2 * n)); % Creates a row vector for the signal x comprising 2 sinusoids
X = fftshift(fft(x, Nzp)); % computes the DFT of DT signal x with zero-padding up to Nzp samples.

w = hamming(N); %creates Hamming Window vector as a column vector
w = w'; % transposes Hamming Window vector so as to obtain a row vector
xw = x .* w; % uses Hamming Window to create windowed signal 
XW = fftshift(fft(xw, Nzp));  % computes the DFT of windowed DT signal xw with zero-padding

freq = ((-Nzp/2):((Nzp/2)-1)) * (Fs/Nzp);  % creates frequency points for plot in Hz
%% Plots 

%Plotting the DFT of x without the Hamming Window
figure 
subplot(2,1,1)
plot(freq(freq>0), 20*log10(abs(X(freq > 0))));
grid on 
xlabel('Frequency (Hz)')
ylabel('|X(\Omega)| (dB)')
title('DFT of 2 sinusoids computed without Hamming Window')

% Plotting the DFT of xw (with Hamming Window used)
subplot(2,1,2)
plot(freq(freq>0), 20*log10(abs(XW(freq>0))));
grid on 
xlabel('Frequency (Hz)')
ylabel('|X(\Omega)| (dB)')
title('DFT of Windowed Input Signal')

end

