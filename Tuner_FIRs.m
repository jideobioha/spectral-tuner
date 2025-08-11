function [b_E, b_A, b_D] = Tuner_FIRs()
%  Tuner_FIRs -     This function will use the Parks-McClellan procedure
%                   to design the filters and compute the FIR filter 
%                   coefficient vectors that are then provided as outputs of the 
%                   function 
%
%   USAGE:   [b_E, b_A, b_D] = Tuner_FIRs();
%   
%   inputs:   No inputs
%
%   outputs:  This function outputs three sets of filter coefficients
%             called b_E, b_A, and b_D. 

Fs = 4410;
rp = 1; rs = 60; % specify passband 
AA = [0 1 0]; %%% specfies that you want a band pass filter (for strings E, A, and D)
dev=[10^(-rs/20) (10^(rp/20)-1)/(10^(rp/20)+1) 10^(-rs/20)]; % parm. needed by design routine for a band pass filter

Omega = 0:0.001:pi; % creates a vector of DT frequency points
freq = Omega * ((Fs)/(2*pi));
%% designing band pass filter for Open Low E string
% % Band Pass Edges
f_spec_E = [59 68 100 109];
[N_E, fo_E, ao_E, w_E] = firpmord(f_spec_E, AA, dev, Fs);
b_E = firpm(N_E, fo_E, ao_E, w_E); % Computes the designed filter coefficients in vector b_E

%% designing band pass filter for Open A string
% % Band Pass Edges
f_spec_A = [86 100 132 141];
[N_A, fo_A, ao_A, w_A] = firpmord(f_spec_A, AA, dev, Fs);
b_A = firpm(N_A, fo_A, ao_A, w_A); % Computes the designed filter coefficients in vector b_A

%% designing band pass filter for Open D string
% % Band Pass Edges
f_spec_D = [120 132 165 174];
[N_D, fo_D, ao_D, w_D] = firpmord(f_spec_D, AA, dev, Fs);
b_D = firpm(N_D, fo_D, ao_D, w_D); % Computes the designed filter coefficients in vector b_A
%% filter plots
figure;
H_E = freqz(b_E, 1, Omega);
H_A = freqz(b_A, 1, Omega);
H_D = freqz(b_D, 1, Omega);
subplot(3,1,1)
plot(freq, 20*log10(abs(H_E)),'b');
xlabel("Frequency (Hz)"); ylabel("|H(\Omega)|_{E} (dB)"); title("Magnitude of Frequency Responses of E band filter")
grid on 
subplot(3,1,2)
plot(freq, 20*log10(abs(H_A)), 'r');
xlabel("Frequency (Hz)"); ylabel("|H(\Omega)|_{A} (dB)"); title("Magnitude of Frequency Responses of A band filter")
grid on 
subplot(3,1,3)
plot(freq, 20*log10(abs(H_D)), 'g');
xlabel("Frequency (Hz)"); ylabel("|H(\Omega)|_{D} (dB)"); title("Magnitude of Frequency Responses of D band filter")
grid on
end