% ================================================================
% compare.m
%
% Reconstructs the ternary sigma-delta stream and compares it with
% the floating-point FDTD reference.
%
% Usage:
%
%   octave compare.m                     % uses probe_fixed.txt
%   octave compare.m probe_fixed.txt
%
% or from the Octave prompt:
%
%   fn = 'probe_fixed.txt'; compare
%
% No toolbox packages are required.  sinc, the Hann window, the
% 'same' convolution and the Pearson correlation are all computed
% inline so that neither the signal nor the statistics package is
% needed.
%
% Expected input columns:
%
%   step   ternary_symbol   fp_Ez   [state_Ez]
%
% ================================================================

1;

% ------------------------------------------------
% Helper: sinc(x) = sin(pi x)/(pi x), sinc(0) = 1
% ------------------------------------------------
function y = sinc_(x)

    y = ones(size(x));

    nz = (x ~= 0);

    y(nz) = sin(pi*x(nz)) ./ (pi*x(nz));

end

% ------------------------------------------------
% Helper: symmetric Hann window, identical to
% numpy.hanning(N)
% ------------------------------------------------
function w = hann_(N)

    if N == 1
        w = 1;
        return;
    end

    w = 0.5 * (1 - cos(2*pi*(0:N-1)'/(N-1)));

end

% ------------------------------------------------
% Helper: 'same' convolution, identical to
% numpy.convolve(x, h, 'same') for odd length(h)
% ------------------------------------------------
function y = conv_same(x, h)

    x = x(:);
    h = h(:);

    N = length(x);
    M = length(h);

    full = conv(x, h);

    start = floor((M-1)/2) + 1;

    y = full(start : start + N - 1);

end

% ------------------------------------------------
% Helper: population standard deviation
% (numpy default, ddof = 0)
% ------------------------------------------------
function s = pstd(x)

    s = std(x(:), 1);

end

% ------------------------------------------------
% Helper: Pearson correlation coefficient
% ------------------------------------------------
function r = pearson(a, b)

    a = a(:) - mean(a(:));
    b = b(:) - mean(b(:));

    denom = sqrt(sum(a.^2) * sum(b.^2));

    if denom == 0
        r = NaN;
    else
        r = sum(a .* b) / denom;
    end

end

% ------------------------------------------------
% Helper: normalized one-sided magnitude spectrum
% ------------------------------------------------
function S = spec_(x)

    x = x(:);

    N = length(x);

    x = (x - mean(x)) .* hann_(N);

    X = abs(fft(x));

    S = X(1 : floor(N/2) + 1);

    S = S / max(S);

end

% ------------------------------------------------
% Helper: report block
% ------------------------------------------------
function report(a, b, label)

    a = a(:);
    b = b(:);

    err = norm(a - b) / norm(b);

    printf('%-34s Pearson = % .6f   NRMSE = %6.3f %%   amp ratio = %.4f\n', ...
           label, ...
           pearson(a, b), ...
           err*100, ...
           pstd(a)/pstd(b));

end

% ================================================================
% Main
% ================================================================

AMP = 20.0;

% ------------------------------------------------
% Input file
% ------------------------------------------------

if ~exist('fn', 'var')

    args = argv();

    if numel(args) >= 1
        fn = args{1};
    else
        fn = 'probe_fixed.txt';
    end

end

if exist(fn, 'file') ~= 2
    error('Cannot find input file: %s', fn);
end

% Skip the one-line text header.
d = dlmread(fn, '', 1, 0);

if columns(d) < 3
    error('Input file must contain at least three columns.');
end

n   = d(:,1);
sym = d(:,2);
fp  = d(:,3);

% ------------------------------------------------
% Reconstruction filter
%
% All FDTD cavity modes lie below
%
%   f = asin(CFL*sqrt(2))/pi ~ 0.045 cycles/step
%
% so a linear-phase FIR with cutoff 0.06 passes the
% entire physical band and rejects the shaped
% quantisation noise.
% ------------------------------------------------

NT = 401;
fc = 0.06;

k = (0:NT-1)' - (NT-1)/2;

h = sinc_(2*fc*k) .* hann_(NT);
h = h / sum(h);

D = floor((NT-1)/2);

rec = AMP * conv_same(sym, h);      % reconstructed ternary field
fpf =       conv_same(fp,  h);      % same filter on the reference

% ------------------------------------------------
% Discard filter edges and startup transient
%
% Python slice [s0 : s1] with s0 = 600 (0-based) and
% s1 = numel(n) - NT  becomes 1-based s0+1 : s1 here.
% ------------------------------------------------

s0 = 600;
s1 = numel(n) - NT;

idx = (s0+1) : s1;

r = rec(idx);
f = fpf(idx);
t = n(idx);

% ------------------------------------------------
% Summary
% ------------------------------------------------

printf('\n');
printf('samples compared        : %d\n', numel(t));

printf('ternary symbol usage    : -1:%.1f%%  0:%.1f%%  +1:%.1f%%\n', ...
       mean(sym == -1)*100, ...
       mean(sym ==  0)*100, ...
       mean(sym ==  1)*100);

printf('FP peak |Ez| at probe   : %.3f  (full scale AMP = %g)\n\n', ...
       max(abs(fp)), AMP);

report(r, f, 'reconstructed vs filtered FP');
report(AMP*sym(idx), fp(idx), 'raw symbols vs raw FP');

% ------------------------------------------------
% Spectra
% ------------------------------------------------

Nt = numel(t);

fr = (0 : floor(Nt/2))' / Nt;       % numpy.fft.rfftfreq(Nt)

St = spec_(r);
Sf = spec_(f);

[~, it] = max(St(2:end));
[~, iff] = max(Sf(2:end));

printf('\n');
printf('dominant freq  ternary  : %.6f cycles/step\n', fr(it+1));
printf('dominant freq  FP       : %.6f cycles/step\n', fr(iff+1));
printf('source frequency        : 0.010000 cycles/step\n');

band = (fr <= 0.06);

cosim = dot(St(band), Sf(band)) / (norm(St(band)) * norm(Sf(band)));

printf('in-band spectral cosine similarity : %.6f\n', cosim);

% ------------------------------------------------
% Where the quantisation noise went
% ------------------------------------------------

Sraw = spec_(sym(idx));

hi = (fr > 0.1);

printf('raw stream: fraction of spectral power above f=0.1 : %.1f %%\n', ...
       sum(Sraw(hi).^2) / sum(Sraw.^2) * 100);

% ------------------------------------------------
% Save the reconstructed traces
% ------------------------------------------------

arr = [t(:), r(:), f(:)];

save('-ascii', 'arr.txt', 'arr');

% ================================================================
% SPECTRAL COMPARISON
%
% Everything above is unchanged. This section adds amplitude-
% corrected spectra, quantitative in-band comparison, and plots.
% ================================================================

CFL = 0.1;

% Highest frequency any cavity mode of this grid can occupy.
% Beyond this the discrete wave operator has no support, so
% anything there is quantisation noise.
f_mode_max = asin(CFL*sqrt(2))/pi;

% ------------------------------------------------
% Amplitude-corrected single-sided spectra
%
% The normalised spectra St and Sf computed above are fine for
% shape comparison but hide the scale. These are in field units.
% ------------------------------------------------

w  = hann_(Nt);
cg = mean(w);                        % Hann coherent gain = 0.5

Xr = fft((r - mean(r)) .* w);
Xf = fft((f - mean(f)) .* w);
Xs = fft((sym(idx) - mean(sym(idx))) .* w);

nb = floor(Nt/2) + 1;

Ar = 2*abs(Xr(1:nb)) / (Nt*cg);      % reconstructed ternary
Af = 2*abs(Xf(1:nb)) / (Nt*cg);      % filtered reference
As = 2*abs(Xs(1:nb)) / (Nt*cg);      % raw symbol stream

Ar(1) = Ar(1)/2;
Af(1) = Af(1)/2;
As(1) = As(1)/2;

% ------------------------------------------------
% Quantitative comparison
% ------------------------------------------------

[pr, ir] = max(Ar(2:end));  ir = ir + 1;
[pf, iff2] = max(Af(2:end)); iff2 = iff2 + 1;

df = 1/Nt;

printf('\n');
printf('=========================================================\n');
printf(' SPECTRAL COMPARISON\n');
printf('=========================================================\n');
printf('\n');
printf('FFT resolution                 : %.6f cycles/step\n', df);
printf('Highest cavity mode frequency  : %.6f cycles/step\n', f_mode_max);
printf('\n');
printf('Peak frequency  ternary        : %.6f cycles/step\n', fr(ir));
printf('Peak frequency  FP             : %.6f cycles/step\n', fr(iff2));
printf('Peak separation                : %.2f FFT bins\n', abs(fr(ir)-fr(iff2))/df);
printf('\n');
printf('Peak amplitude  ternary        : %.4f\n', pr);
printf('Peak amplitude  FP             : %.4f\n', pf);
printf('Amplitude ratio (ternary/FP)   : %.4f\n', pr/pf);
printf('\n');

% Pearson correlation of the two spectra over the physical band
% only. Comparing across the full band is dominated by the shaped
% noise and tells you nothing about the physics.
inband = (fr <= f_mode_max);

printf('In-band (f <= %.4f) comparison of normalised spectra:\n', f_mode_max);
printf('  Pearson                      : %.6f\n', pearson(St(inband), Sf(inband)));
printf('  cosine similarity            : %.6f\n', ...
       dot(St(inband), Sf(inband)) / (norm(St(inband))*norm(Sf(inband))));
printf('  relative shape error         : %.6f\n', ...
       norm(St(inband)-Sf(inband)) / norm(Sf(inband)));
printf('\n');

% In-band signal-to-quantisation-noise ratio: power within two
% bins of the drive peak against everything else below the mode
% band edge.
sig  = abs(fr - fr(iff2)) <= 2*df;
nois = inband & ~sig;

snr_t = 10*log10(sum(Ar(sig).^2) / sum(Ar(nois).^2));
snr_f = 10*log10(sum(Af(sig).^2) / sum(Af(nois).^2));

printf('In-band SNR  ternary           : %8.2f dB\n', snr_t);
printf('In-band SNR  FP                : %8.2f dB\n', snr_f);
printf('Degradation from quantisation  : %8.2f dB\n', snr_f - snr_t);
printf('\n');

% Where the raw stream keeps its energy.
printf('Raw stream power below mode band edge : %.2f %%\n', ...
       100*sum(As(inband).^2)/sum(As.^2));
printf('Raw stream power above f = 0.1        : %.2f %%\n', ...
       100*sum(As(fr>0.1).^2)/sum(As.^2));
printf('\n');

% ------------------------------------------------
% Figure 1: time domain
% ------------------------------------------------

f1 = figure(1);
clf;
set(f1, 'Position', [60 60 900 400]);

plot(t, f, 'k-',  'LineWidth', 1.6);
hold on;
plot(t, r, 'r-',  'LineWidth', 1.0);
hold off;

grid on;
xlabel('Time step');
ylabel('E_z at (32,32)');
title('Reconstructed ternary field against floating-point reference');
legend('floating point (filtered)', 'ternary (reconstructed)', ...
       'Location', 'northwest');

print(f1, 'compare_time.png', '-dpng', '-r120');

% ------------------------------------------------
% Figure 2: in-band spectra, linear
% ------------------------------------------------

f2 = figure(2);
clf;
set(f2, 'Position', [80 80 900 400]);

plot(fr(inband), Af(inband), 'k-', 'LineWidth', 1.6);
hold on;
plot(fr(inband), Ar(inband), 'r-', 'LineWidth', 1.0);
hold off;

grid on;
xlabel('Frequency (cycles / simulation step)');
ylabel('Amplitude (field units)');
title('Spectra over the physical band');
legend('floating point', 'ternary (reconstructed)', 'Location', 'northeast');
xlim([0 f_mode_max]);

print(f2, 'compare_spectrum_inband.png', '-dpng', '-r120');

% ------------------------------------------------
% Figure 3: full band, log scale
%
% This is the plot that shows the method working: the raw stream
% carries most of its power near Nyquist, far above the highest
% frequency the grid can support, which is why a low-pass
% reconstruction recovers the field at all.
% ------------------------------------------------

f3 = figure(3);
clf;
set(f3, 'Position', [100 100 900 440]);

semilogy(fr, As + eps, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);
hold on;
semilogy(fr, Af + eps, 'k-', 'LineWidth', 1.6);
semilogy(fr, Ar + eps, 'r-', 'LineWidth', 1.0);

yl = ylim();
plot([f_mode_max f_mode_max], yl, 'b--', 'LineWidth', 1.2);
plot([fc fc], yl, 'g--', 'LineWidth', 1.2);
ylim(yl);
hold off;

grid on;
xlabel('Frequency (cycles / simulation step)');
ylabel('Amplitude (field units)');
%title('Full-band spectra: noise shaping pushes error above the mode band');
legend('raw ternary symbols', ...
       'floating point', ...
       'ternary (reconstructed)', ...
       'highest cavity mode', ...
       'reconstruction cutoff', ...
       'Location', 'southeast');
xlim([0 0.5]);

print(f3, 'compare_spectrum_full.png', '-dpng', '-r120');

printf('Written: compare_time.png, compare_spectrum_inband.png, compare_spectrum_full.png\n');
