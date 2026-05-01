%this code is to simulate the gfsk modulation for bt values BT=0.3,0.5,0.7,1 with AWGN and Rayleigh Channels.
%It outputs two graphs, BER(Bit Error Rate) and PSD(Power Spectral Density).
clear; close all; clc;

%--- 1. Simulation Parameters ---
Rb = 1e6;             % Bit rate (1 Mbps)
sps = 8;              % Samples per symbol 
fs = Rb * sps;        % Sampling frequency
h = 0.5;              % Modulation index (nRF24L01+ standard)
filterSpan = 6;       % Gaussian filter length in symbols

%--- Simulation Controls ---
snrDbVec = 0:2:20;               % SNR (Eb/No) range in dB
numBits = 1.5e6;                 % 1.5 Million bits for smooth statistical curves
btValues = [0.3, 0.5, 0.7, 1.0]; % BT Products to test

colors = {'r', 'b', 'g', 'm'};
lines_awgn = {'-o', '-s', '-^', '-d'};
lines_rayleigh = {'--o', '--s', '--^', '--d'};

ber_awgn = zeros(length(btValues), length(snrDbVec));
ber_rayleigh = zeros(length(btValues), length(snrDbVec));
psd_store = cell(length(btValues), 1);
freqVec_store = cell(length(btValues), 1);

bits = randi([0 1], numBits, 1); 

fprintf('Starting Corrected GFSK simulation (Slow Fading Channel)...\n');
mainTic = tic;

for btIdx = 1:length(btValues)
    BT = btValues(btIdx);
    fprintf('  Processing BT = %.1f...\n', BT);
    
    %----------------------------------------------------------------------
    % 2. Perfect Continuous Phase Modulator
    %----------------------------------------------------------------------
    % A. Gaussian Filter
    alpha = sqrt(log(2)/2) / BT;
    t_norm = -filterSpan/2 : 1/sps : filterSpan/2; 
    h_gauss = (sqrt(pi)/alpha) * exp(-(pi * t_norm / alpha).^2);
    h_gauss = h_gauss / sum(h_gauss);

    % B. Pulse Shaping & Filter (ISI Generation)
    nrz_rect = repelem(2*bits - 1, sps); 
    freq_dev = conv(nrz_rect, h_gauss, 'same'); 
    
    % C. Integrate to Phase
    phase = cumsum(freq_dev) * (h * pi / sps);
    txSig = exp(1j * phase);

    %----------------------------------------------------------------------
    % 3. Smooth PSD via Welch's Method
    %----------------------------------------------------------------------
    windowSize = sps * 512;  
    overlap = windowSize / 2;
    nfft = sps * 1024;
    
    [pxx, f] = pwelch(txSig, hamming(windowSize), overlap, nfft, fs, 'centered');
    psd_store{btIdx} = 10*log10(pxx / max(pxx));
    freqVec_store{btIdx} = f / Rb;
    

    ma_len = 40 * sps; 
    h_fade_raw = (randn(length(txSig), 1) + 1j*randn(length(txSig), 1)) / sqrt(2);
    h_fade_smooth = filter(ones(ma_len,1)/ma_len, 1, h_fade_raw);
    
    % Normalize average channel power to 0 dB (1.0 linear)
    h_fade_smooth = h_fade_smooth / sqrt(mean(abs(h_fade_smooth).^2)); 
    
    txSig_faded = txSig .* h_fade_smooth;

    for snrIdx = 1:length(snrDbVec)
        EbNo_dB = snrDbVec(snrIdx);
        
        snr_db_baseband = EbNo_dB - 10*log10(sps);
        noise_var = 10^(-snr_db_baseband/10);
        noise = sqrt(noise_var/2) * (randn(size(txSig)) + 1j*randn(size(txSig)));
        
        rxSig_awgn = txSig + noise;
        rxSig_rayleigh = txSig_faded + noise;
        
        %--- 5. 1-Symbol Delay Differential Demodulator ---
        % AWGN Demod
        rx_delayed_awgn = circshift(rxSig_awgn, sps);
        rx_delayed_awgn(1:sps) = 0; 
        phase_diff_awgn = angle(rxSig_awgn .* conj(rx_delayed_awgn));
        
        % Rayleigh Demod (Now works perfectly because channel phase is stable across 1 symbol)
        rx_delayed_rayleigh = circshift(rxSig_rayleigh, sps);
        rx_delayed_rayleigh(1:sps) = 0; 
        phase_diff_rayleigh = angle(rxSig_rayleigh .* conj(rx_delayed_rayleigh));
        
        sample_indices = sps : sps : (numBits*sps);
        
        demod_bits_awgn = phase_diff_awgn(sample_indices) > 0;
        demod_bits_rayleigh = phase_diff_rayleigh(sample_indices) > 0;
        
        ber_awgn(btIdx, snrIdx) = sum(bits ~= demod_bits_awgn(:)) / numBits;
        
        err_rayleigh = sum(bits ~= demod_bits_rayleigh(:)) / numBits;
        ber_rayleigh(btIdx, snrIdx) = max(err_rayleigh, 1/(numBits*10)); 
    end
end
fprintf('Simulation complete. Total time: %.1f seconds.\n', toc(mainTic));

%==========================================================================
% 6. Plotting Formatting
%==========================================================================

% --- Graph 1: BER ---
figure('Position', [100, 100, 900, 550]);
hold on; grid on; box on;
title('Graph 1: BER Performance Trade-off (GFSK in AWGN & Slow Fading)');
xlabel('SNR (Eb/No in dB)');
ylabel('Bit Error Rate (BER)');
set(gca, 'YScale', 'log'); 
axis([min(snrDbVec) max(snrDbVec) 1e-5 1]); 

legend_str = {}; plot_handles = [];
for btIdx = 1:length(btValues)
    p1 = plot(snrDbVec, ber_awgn(btIdx, :), lines_awgn{btIdx}, 'Color', colors{btIdx}, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', colors{btIdx}, 'MarkerEdgeColor', 'k');
    p2 = plot(snrDbVec, ber_rayleigh(btIdx, :), lines_rayleigh{btIdx}, 'Color', colors{btIdx}, 'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', colors{btIdx});
    legend_str = [legend_str, sprintf('BT=%.1f (AWGN)', btValues(btIdx)), sprintf('BT=%.1f (Rayleigh)', btValues(btIdx))];
    plot_handles = [plot_handles, p1, p2];
end
legend(plot_handles, legend_str, 'Location', 'southwest', 'NumColumns', 2);
set(gcf, 'Color', 'w');

% --- Graph 2: PSD ---
figure('Position', [150, 150, 900, 550]);
hold on; grid on; box on;
title('Graph 2: GFSK Power Spectral Density (Smoothed)');
xlabel('Normalized Frequency (f - f_c) / R_b');
ylabel('Power Spectral Density (dB/Hz)');
axis([-2 2 -80 5]); 

psd_styles = {'-', '--', '-.', ':'};
for btIdx = 1:length(btValues)
    plot(freqVec_store{btIdx}, psd_store{btIdx}, psd_styles{btIdx}, 'Color', colors{btIdx}, 'LineWidth', 2.0);
end
legend({'BT=0.3', 'BT=0.5', 'BT=0.7', 'BT=1.0'}, 'Location', 'northeast');
set(gcf, 'Color', 'w');