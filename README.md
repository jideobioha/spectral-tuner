# spectral-tuner
EECE 301 - Signals and Systems - DFT-based guitar tuner

DFT-Based Guitar Tuner
A simple guitar tuner built using the Discrete Fourier Transform (DFT) to estimate a string’s fundamental frequency with ±5 cents accuracy.

What It Does
Analyzes audio signals from real or synthetic guitar strings

Detects fundamental frequency and compares to standard tuning

Identifies closest note for open E, A, and D strings

How It Works
Signal Analysis – Explored effects of sampling rate, zero-padding, and windowing (Hamming) on DFT accuracy.

Harmonic Detection – Verified Fourier Series predictions using guitar recordings.

Synthetic Notes – Generated clean test tones using Fourier Series modeling.

Tuner Design – Used FIR bandpass filters (Parks–McClellan) to isolate note ranges and determine which string is being played.

Key Features
Works with real and synthetic signals

Frequency estimation within ±5 cents

Harmonic and spectral analysis tools

Limitations
Only detects open E, A, D strings

Lower limit detection at 68 Hz is less reliable

 