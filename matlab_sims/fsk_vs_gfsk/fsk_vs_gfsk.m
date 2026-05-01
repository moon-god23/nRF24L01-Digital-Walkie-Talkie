%this code is to simulate the fsk and gfsk modulation for bt values BT=0.5 and BT=1
%it outputs two graphs BER(Bit Error Rate) and PSD(Power Spectral Density)

clear; close all; clc;

%--- 1. Simulation Parameters ---
Rb = 1e6;             % Bit rate (1 Mbps)
sps = 8;              % Samples per symbol 
fs = Rb * sps;        % Sampling frequency
h = 0.5;              % Modulation index
filterSpan = 6;       % Gaussian filter length in symbols

%--- Simulation Controls ---
snrDbVec = 0:2:20;               % SNR (Eb/No) range in dB
numBits = 1.5e6;                 % 1.5 Million bits for smooth curves
btValues = [0.5, 1.0];           % GFSK BT Products to test

%--- Storage Arrays ---
ber_fsk_awgn  = zeros(1, length(snrDbVec));
ber_fsk_ray   = zeros(1, length(snrDbVec));
ber_gfsk_awgn = zeros(length(btValues), length(snrDbVec));
ber_gfsk_ray  = zeros(length(btValues), length(snrDbVec));

psd_store = cell(3, 1);       
freqVec_store = cell(3, 1);

fprintf('Generating bits...\n');
bits = randi([0 1], 1, numBits); 
nrz_rect = repelem(2*bits - 1, sps); 

%==========================================================================
% 2. Standard FSK & GFSK Generation
%==========================================================================
fprintf('Processing Modulations...\n');

% --- Standard FSK ---
phase_fsk = cumsum(nrz_rect) * (h * pi / sps);
txSig_fsk = exp(1j * phase_fsk);

windowSize = sps * 512; overlap = windowSize / 2; nfft = sps * 1024;
[pxx, f] = pwelch(txSig_fsk, hamming(windowSize), overlap, nfft, fs, 'centered');
psd_store{1} = 10*log10(pxx / max(pxx));
freqVec_store{1} = f / Rb;

% --- GFSK (BT = 0.5 and 1.0) ---
txSig_gfsk = zeros(length(btValues), length(txSig_fsk));

for btIdx = 1:length(btValues)
    BT = btValues(btIdx);
    
    alpha = sqrt(log(2)/2) / BT;
    t_norm = -filterSpan/2 : 1/sps : filterSpan/2; 
    h_gauss = (sqrt(pi)/alpha) * exp(-(pi * t_norm / alpha).^2);
    h_gauss = h_gauss / sum(h_gauss);
    
    freq_dev = conv(nrz_rect, h_gauss, 'same'); 
    phase_gfsk = cumsum(freq_dev) * (h * pi / sps);
    txSig_gfsk(btIdx, :) = exp(1j * phase_gfsk);
    
    [pxx, f] = pwelch(txSig_gfsk(btIdx, :), hamming(windowSize), overlap, nfft, fs, 'centered');
    psd_store{btIdx + 1} = 10*log10(pxx / max(pxx));
    freqVec_store{btIdx + 1} = f / Rb;
end

%==========================================================================
% 3. Rayleigh Fading Channel Generation (Slow Fading)
%==========================================================================
fprintf('Generating Rayleigh Fading Profile...\n');
ma_len = 40 * sps; 
h_fade_raw = (randn(1, length(txSig_fsk)) + 1j*randn(1, length(txSig_fsk))) / sqrt(2);
h_fade_smooth = filter(ones(1, ma_len)/ma_len, 1, h_fade_raw);

h_fade_smooth = h_fade_smooth / sqrt(mean(abs(h_fade_smooth).^2)); 

txSig_fsk_faded = txSig_fsk .* h_fade_smooth;
txSig_gfsk_faded = zeros(size(txSig_gfsk));
for btIdx = 1:length(btValues)
    txSig_gfsk_faded(btIdx, :) = txSig_gfsk(btIdx, :) .* h_fade_smooth;
end

%==========================================================================
% 4. Channel Noise & Demodulation
%==========================================================================
fprintf('Running AWGN & Rayleigh Simulation (This may take a minute)...\n');
sample_indices = sps : sps : (numBits*sps);

for snrIdx = 1:length(snrDbVec)
    EbNo_dB = snrDbVec(snrIdx);
    snr_db_baseband = EbNo_dB - 10*log10(sps);
    noise_var = 10^(-snr_db_baseband/10);
    
    noise = sqrt(noise_var/2) * (randn(1, length(txSig_fsk)) + 1j*randn(1, length(txSig_fsk)));
    
    % --- FSK Demod ---
    % AWGN (FIXED: Added sample_indices downsampling)
    rx_fsk_a = txSig_fsk + noise;
    rx_del_fsk_a = circshift(rx_fsk_a, sps); rx_del_fsk_a(1:sps) = 0; 
    demod_fsk_a = angle(rx_fsk_a .* conj(rx_del_fsk_a));
    ber_fsk_awgn(snrIdx) = sum(bits ~= (demod_fsk_a(sample_indices) > 0)) / numBits;
    
    % Rayleigh (FIXED: Added sample_indices downsampling)
    rx_fsk_r = txSig_fsk_faded + noise;
    rx_del_fsk_r = circshift(rx_fsk_r, sps); rx_del_fsk_r(1:sps) = 0; 
    demod_fsk_r = angle(rx_fsk_r .* conj(rx_del_fsk_r));
    ber_fsk_ray(snrIdx) = sum(bits ~= (demod_fsk_r(sample_indices) > 0)) / numBits;
    
    % --- GFSK Demod ---
    for btIdx = 1:length(btValues)
        % AWGN
        rx_gfsk_a = txSig_gfsk(btIdx, :) + noise;
        rx_del_gfsk_a = circshift(rx_gfsk_a, sps); rx_del_gfsk_a(1:sps) = 0; 
        demod_a = angle(rx_gfsk_a .* conj(rx_del_gfsk_a));
        ber_gfsk_awgn(btIdx, snrIdx) = sum(bits ~= (demod_a(sample_indices) > 0)) / numBits;
        
        % Rayleigh
        rx_gfsk_r = txSig_gfsk_faded(btIdx, :) + noise;
        rx_del_gfsk_r = circshift(rx_gfsk_r, sps); rx_del_gfsk_r(1:sps) = 0; 
        demod_r = angle(rx_gfsk_r .* conj(rx_del_gfsk_r));
        ber_gfsk_ray(btIdx, snrIdx) = sum(bits ~= (demod_r(sample_indices) > 0)) / numBits;
    end
end
fprintf('Simulation complete.\n');

%==========================================================================
% 5. Plotting
%==========================================================================

% --- Graph 1: BER ---
figure('Position', [100, 100, 900, 600]);
hold on; grid on; box on;
title('Graph 1: BER Performance (FSK vs GFSK) in AWGN & Rayleigh Channels');
xlabel('SNR (Eb/No in dB)');
ylabel('Bit Error Rate (BER)');
set(gca, 'YScale', 'log'); 
axis([min(snrDbVec) max(snrDbVec) 1e-5 1]); 

% AWGN Lines (Solid)
p1 = plot(snrDbVec, ber_fsk_awgn, '-ko', 'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', 'k');
p2 = plot(snrDbVec, ber_gfsk_awgn(1, :), '-b^', 'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', 'b');
p3 = plot(snrDbVec, ber_gfsk_awgn(2, :), '-rs', 'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', 'r');

% Rayleigh Lines (Dashed)
p4 = plot(snrDbVec, ber_fsk_ray, '--ko', 'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', 'w');
p5 = plot(snrDbVec, ber_gfsk_ray(1, :), '--b^', 'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', 'w');
p6 = plot(snrDbVec, ber_gfsk_ray(2, :), '--rs', 'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', 'w');

legend([p1, p2, p3, p4, p5, p6], ...
    {'FSK (AWGN)', 'GFSK BT=0.5 (AWGN)', 'GFSK BT=1.0 (AWGN)', ...
     'FSK (Rayleigh)', 'GFSK BT=0.5 (Rayleigh)', 'GFSK BT=1.0 (Rayleigh)'}, ...
    'Location', 'southwest', 'FontSize', 10, 'NumColumns', 2);
set(gcf, 'Color', 'w');

% --- Graph 2: PSD ---
figure('Position', [150, 150, 800, 500]);
hold on; grid on; box on;
title('Graph 2: Power Spectral Density (Spectral Splatter Comparison)');
xlabel('Normalized Frequency (f - f_c) / R_b');
ylabel('Power Spectral Density (dB/Hz)');
axis([-2 2 -80 5]); 

plot(freqVec_store{1}, psd_store{1}, '-k', 'LineWidth', 1.5);
plot(freqVec_store{2}, psd_store{2}, '-b', 'LineWidth', 2.0);
plot(freqVec_store{3}, psd_store{3}, '-r', 'LineWidth', 2.0);

legend({'Standard FSK (High Splatter)', 'GFSK BT = 0.5 (Optimal)', 'GFSK BT = 1.0'}, 'Location', 'northeast', 'FontSize', 11);
set(gcf, 'Color', 'w');