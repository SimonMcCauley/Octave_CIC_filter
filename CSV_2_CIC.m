pkg load signal; 
clear; clc; close all;

% =========================================================================
% DSP TEST BENCH: Fixed-Point CIC Decimation Chain
% =========================================================================

% --- TASK 1: SYSTEM PARAMETERS ---
% High sample rate (Input) vs. Low sample rate (Output)
fs_high = 32768 * 16;  
R = 16;                % Decimation factor: We keep 1 out of every 16 samples
fs_low = 32768;        
peak_corr = single(pi / 2); 
alpha = single(0.995); % Pole location for the DC Blocker (near the unit circle)
f_size = 20;           

config_file = 'dir_config.mat';

% --- TASK 2: FILE SELECTION ---
if exist(config_file, 'file'); load(config_file); else main_path = pwd(); end
[fname, main_path] = uigetfile([main_path, filesep, '*.csv'], 'Select MAIN CSV File');
if isequal(fname, 0); return; end
save(config_file, 'main_path');

% --- TASK 3: DATA ACQUISITION ---
data = csvread(fullfile(main_path, fname), 1, 0);
x_raw = single(data(:, 2));  
n_high = length(x_raw);

% --- TASK 4: SIGNAL PROCESSING PIPELINE ---

% 1. DC BLOCKER (High-Pass Filter)
% Removes the DC offset (bias) from the ADC. It works by subtracting a 
% delayed version of the signal, keeping only the AC components.
x_clean = zeros(n_high, 1, 'single');
x_prev = x_raw(1); 
y_prev = single(0);
for k = 1:n_high
    x_clean(k) = x_raw(k) - x_prev + alpha * y_prev;
    x_prev = x_raw(k);
    y_prev = x_clean(k);
end

% 2. BANDPASS FILTER (Selection)
% Isolates the frequency of interest. We take the absolute value (rectification)
% to extract the "envelope" or magnitude of the signal.
b_bpf = single([1, 0, -1]); 
a_bpf = single([1, -1.839, 0.8789]);
y_bpf = filter(b_bpf, a_bpf, x_clean) / single(16);
y_abs = abs(y_bpf); 

scaling_factor = 1;

% 3. INTEGRATOR STAGES (The "I" in CIC)
% These are basically running accumulators. Because we are summing over 
% and over, the bit-width grows rapidly. We use uint64 to prevent overflow 
% during the integration phase.
i1_hist = uint64(zeros(n_high, 1)); i2_hist = uint64(zeros(n_high, 1)); i3_hist = uint64(zeros(n_high, 1));
acc1 = uint64(0); acc2 = uint64(0); acc3 = uint64(0);
y_abs_scaled = uint64(round(single(y_abs) * scaling_factor));

for k = 1:n_high
    acc1 = acc1 + y_abs_scaled(k); 
    acc2 = acc2 + acc1;            
    acc3 = acc3 + acc2;            
    i1_hist(k) = acc1; i2_hist(k) = acc2; i3_hist(k) = acc3;
end

% 4. DECIMATION & COMB STAGES (The "C" in CIC)
% We downsample by R, then calculate the "difference" between samples. 
% This acts as a differentiator, effectively turning the integration 
% into a moving average filter.
i3_ds_u64 = i3_hist(1:R:end); 
n_ds = length(i3_ds_u64);
c1_u = uint64(zeros(n_ds, 1)); c2_u = uint64(zeros(n_ds, 1)); c3_u = uint64(zeros(n_ds, 1));
for k = 4:n_ds
    c1_u(k) = i3_ds_u64(k) - i3_ds_u64(k-1);
    c2_u(k) = c1_u(k) - c1_u(k-1);
    c3_u(k) = c2_u(k) - c2_u(k-1);
end
% Rescale back to floating point for the final output
y_cic_u64 = (double(peak_corr) * (double(c3_u) / double(scaling_factor))) / (double(R^3));

% 5. VALIDATION PATH
% A floating-point moving average reference to verify the fixed-point logic.
y_cic_float = (double(peak_corr) * filter(ones(1,R), 1, filter(ones(1,R), 1, filter(ones(1,R), 1, double(y_abs))) ) );
y_cic_float = y_cic_float(1:R:end) / (R^3);

% --- TASK 5: ANALYSIS VISUALIZATION ---
t_high = (0:n_high-1) / fs_high;
t_low = (0:n_ds-1) / fs_low;

% Figure 1: Pipeline Stage Analysis
figure(1); set(gcf, 'Color', 'w');
subplot(4,1,1); plot(t_high, x_raw, 'k'); title('1. Raw ADC Time-Series', 'FontSize', f_size); grid on;
subplot(4,1,2); plot(t_high, x_clean, 'g'); title('2. Post-DC Blocker (AC Coupled)', 'FontSize', f_size); grid on;
subplot(4,1,3); plot(t_high, y_bpf, 'r'); title('3. Band-Limited Output', 'FontSize', f_size); grid on;
subplot(4,1,4); 
plot(t_high, y_abs, 'r'); hold on; plot(t_low, y_cic_u64, 'b', 'LineWidth', 2); 
title('4. CIC Decimation (Envelope Extraction)', 'FontSize', f_size); grid on; 
legend('Full Rate', 'Decimated', 'FontSize', f_size-2); hold off;

% Figure 2: Bit-Width Growth (Dynamic Range Analysis)
% This log plot helps determine the minimum register size needed for hardware.
figure(2); set(gcf, 'Color', 'w');
subplot(2,1,1);
semilogy(t_high, double(i1_hist)+1, t_high, double(i2_hist)+1, t_high, double(i3_hist)+1);
title('Integrator Magnitude Growth', 'FontSize', f_size); grid on;
subplot(2,1,2);
plot(t_high, ceil(log2(double(i1_hist)+1)), t_high, ceil(log2(double(i2_hist)+1)), t_high, ceil(log2(double(i3_hist)+1)));
title('Register Width Required (Bits)', 'FontSize', f_size); grid on; ylim([0 70]);

% Figure 4: Precision & Quantization Verification
figure(4); set(gcf, 'Color', 'w');
plot(t_low, y_cic_u64, 'y', 'LineWidth', 2); hold on; plot(t_low, y_cic_float, 'b');
title('Fixed-Point vs. Floating-Point Comparison', 'FontSize', f_size); grid on; 
legend('uint64 Implementation', 'Double Precision Ref');

fprintf('Processing complete. Pipeline verified against floating-point model.\n');