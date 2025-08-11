function fo_est = DFT_Tuner(x, b_E, b_A, b_D)
%  DFT_Tuner  -  The DFT_Tuner function gives an estimate of the fundamental
%                frequency of a guitar signal by applying 3 bandpass
%                filters, (E band, A band and D band filters). It then
%                computes the power of the signal after passing through
%                each of these filters to identify which band is the
%                dominant one. This function then uses DFT analysis to
%                identify the peak frequency (the fundamental frequency of
%                the guitar signal) and passes this frequency as the output. 
%                If the signal's fundamental frequency is outside the frequency
%                bands of interest, the function outputs NaN
%
%   USAGE:   fo_est = DFT_Tuner(x, b_E, b_A, b_D);
%   
%   inputs:   x - The input audio signal (for a guitar)
%             b_E - bandpass filter for detecting fundamental frequency in
%             the E band
%
%             b_A - bandpass filter for detecting fundamental frequency in
%             the A band
%
%             b_D - bandpass filter for detecting fundamental frequency in
%             the D band
%
%   outputs:  This function outputs the estimated frequency, fo_est of the input 
%             signal passed in. In addition, this function creates a single figure with 
%             two subplots: one showing the full DFT (in dB) over all the 
%             positive frequencies and one showing the DFT (in dB) only over the 
%             range of frequencies of the extracted section used while searching 
%             for the peak


Fs = 4410; %sets the sampling frequency
Nzp = 32768; % amount of zero-padding to be used when computing DFT
freq = ((-Nzp/2):((Nzp/2)-1)) * (Fs/Nzp);  % creates frequency points for plot in Hz
fo_est = 0; % initialize estimated frequency

%% Creating Filters to detect if fundamental frequency of filter is below or above frequency window
rp = 1; rs = 60; % specify passband 

% Creating Low Pass filter with a cut off frequency of 68 Hz to detect if
% frequency is below our window of interest
AA_lpf = [1 0]; %%% specifies that you want an LPF
dev_lpf=[(10^(rp/20)-1)/(10^(rp/20)+1) 10^(-rs/20)]; % parm. needed by design routine for LPF
f_spec_lpf = [68 77]; % specify stopband and passband edges
[N_lpf, fo_lpf, ao_lpf, w_lpf] = firpmord(f_spec_lpf, AA_lpf, dev_lpf, Fs);
b_lpf = firpm(N_lpf, fo_lpf, ao_lpf, w_lpf);


% Creating High Pass filter with a cut off frequency of 165 Hz to detect if
% frequency is above our window of interest
AA_hpf = [0 1]; %%% specfies that you want an LPF
dev_hpf=[10^(-rs/20) (10^(rp/20)-1)/(10^(rp/20)+1)]; % parm. needed by design routine for LPF
f_spec_hpf = [165 175]; % specify stopband and passband edges
[N_hpf, fo_hpf, ao_hpf, w_hpf] = firpmord(f_spec_hpf, AA_hpf, dev_hpf, Fs);
b_hpf = firpm(N_hpf, fo_hpf, ao_hpf, w_hpf);

%% Passing input signal through filters

% Passing input through E-band filter
y_E = filter(b_E, 1, x);

% Passing input through A-band filter
y_A = filter(b_A, 1, x);

% Passing input through D-band filter
y_D = filter(b_D, 1, x);

% Passing input through Low Pass filter
y_lpf = filter(b_lpf, 1, x);

% Passing input through High Pass filter
y_hpf = filter(b_hpf, 1, x);

%% determining which band input signal's fundamental frequency lies in 

% initializing variables that will store the power of each filtered signal
% and the max power of all three signals 
E_filter_power = 0;
A_filter_power = 0;
D_filter_power = 0;
lpf_power = 0;
hpf_power = 0;
max_power_signal = 0;


%N:B y_E, y_A, y_D each have the same length

% computing the sum of the squares in each of the filtered signals
for k = 1:length(y_E)
    E_filter_power = E_filter_power + y_E(k).^2;
    A_filter_power = A_filter_power + y_A(k).^2;
    D_filter_power = D_filter_power + y_D(k).^2;
    lpf_power = lpf_power + y_lpf(k).^2;
    hpf_power = hpf_power + y_hpf(k).^2;
end

% dividing sum of squares by length of filtered signals to obtain power
E_filter_power = E_filter_power / length(y_E);
A_filter_power = A_filter_power / length(y_A);
D_filter_power = D_filter_power / length(y_D);
lpf_power = lpf_power / length(y_lpf);
hpf_power = hpf_power / length(y_hpf);

% determining which frequency band has the most power
max_power_signal = max([E_filter_power, A_filter_power, D_filter_power, lpf_power, hpf_power]);

%% V-C: DFT Processing

% create hamming vector
w = hamming(length(x));
w = w.';

% apply window to input audio vector
xw = x .* w;

%compute DFT 
XW = fftshift(fft(xw, Nzp));

%compute the DFT of filtered signals
Y_E = fftshift(fft(y_E, Nzp));
Y_A = fftshift(fft(y_A, Nzp));
Y_D = fftshift(fft(y_D, Nzp));

% finding the band that the input signal's frequency lies in 
% extracts that band's portion of the DFT and reassign it to XW
%also extracting the frequency range of that band with some allowance
if max_power_signal == E_filter_power
    lower_bound = 65;
    upper_bound = 103;
    XW = XW(freq >= lower_bound & freq <= upper_bound);
    band_freq = freq(freq >= lower_bound & freq <= upper_bound);
    disp(['The fundamental frequency of the input audio signal' ...
    ' is within the E band']);
elseif max_power_signal == A_filter_power
    lower_bound = 97;
    upper_bound = 135;
    XW = XW( freq >= lower_bound & freq <= upper_bound);
    band_freq = freq(freq >= lower_bound & freq <= upper_bound);
    disp(['The fundamental frequency of the input audio signal' ...
    ' is within the A band']);
elseif max_power_signal == D_filter_power
    lower_bound = 129;
    upper_bound = 168;
    XW = XW(freq >= lower_bound & freq <= upper_bound);
    band_freq = freq(freq >= lower_bound & freq <= upper_bound);
    disp(['The fundamental frequency of the input audio signal' ...
    ' is within the D band']);
end

% checking for if the input signal's fundamental frequency range lies out
% of bound
if max_power_signal == lpf_power || max_power_signal == hpf_power
    % set the estimated frequency to invalid if singal's fundamental frequency
    % is outside our window of interest
    fo_est = NaN; 

    % indicating which side invalid frequency lies on
    if max_power_signal == lpf_power
        disp(['The fundamental frequency of the input audio signal' ...
            ' is below 68Hz']);
    else
        disp(['The fundamental frequency of the input audio signal' ...
            ' is above 165 Hz']);
    end

else
    % finding the peak of the DFT passed through
    max_dft = max(20*log10(abs(XW)));
  
    % using index to get frequency that corresponds to max dft
    fo_est = band_freq(20*log10(abs(XW)) == max_dft);

    % displaying the fundamental frequency of the signal
    disp(['The fundamental frequency of the audio signal is approximately ' num2str(fo_est) '.'])
end


%% V-D Plots - only plotting for fundamental frequencies within range of interest

if ~isnan(fo_est)

    % Plotting DFT of signal over positive frequencies
    XW_full = fftshift(fft(xw, Nzp));
    figure;
    subplot(2,1,1)
    plot(freq(freq > 0), 20*log10(abs(XW_full(freq > 0))));
    xlabel('Frequency (Hz)'); ylabel('|X(\Omega)| (dB)');
    title('Full DFT of Audio Signal');
    grid on
    
    % Plotting the 'passed' portion of input signal's DFT
    subplot(2,1,2)
    plot(freq(freq >= lower_bound & freq <= upper_bound), 20*log10(abs(XW)));
    xlabel('Frequency (Hz)'); ylabel('|X(\Omega)| (dB)');
    title("DFT of 'passed' portion of Audio Signal");
    grid on 
end

%% Plotting low pass and high pass filters for out-of-range notes
% plotting Low Pass Filter
Omega = 0:0.001:pi; 
H_lpf = freqz(b_lpf, 1, Omega);
ff = Omega * Fs/(2*pi);
figure; 
subplot(2,1,1)
plot(ff, 20*log10(abs(H_lpf)), "m");
xlabel('Frequency (Hz)'); ylabel('|X(\Omega)_{LPF}|   (dB)');
title("Low Pass Filter with a Cutoff Frequency of 68 Hz");

%plotting High Pass Filter
H_hpf = freqz(b_hpf, 1, Omega);
subplot(2,1,2)
plot(ff, 20*log10(abs(H_hpf)), "r");
xlabel('Frequency (Hz)'); ylabel('|H(\Omega)_{HPF}|   (dB)');
title("High Pass Filter with a Cutoff Frequency of 165 Hz");

% plotting DFT of filtered signals
figure
subplot(3,1,1)
plot(freq, 20*log10(abs(Y_E)), 'k');
xlabel('Frequency (Hz)'); ylabel('|X(\Omega)_{LPF}|   (dB)');
title('DFT of input signal after being filtered through E-band filter');
axis([0 1000 -80 50])

subplot(3,1,2)
plot(freq, 20*log10(abs(Y_A)));
xlabel('Frequency (Hz)'); ylabel('|X(\Omega)_{LPF}|   (dB)');
title('DFT of input signal after being filtered through A-band filter');
axis([0 1000 -80 50])

subplot(3,1,3)
plot(freq, 20*log10(abs(Y_D)), 'm');
xlabel('Frequency (Hz)'); ylabel('|X(\Omega)_{LPF}|   (dB)');
title('DFT of input signal after being filtered through D-band filter');
axis([0 1000 -80 50])
end