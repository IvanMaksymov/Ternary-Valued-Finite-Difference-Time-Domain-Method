% ================================================================
% snapshot_overlay.m
%
% Compares ternary and floating-point Ez snapshots.
%
% Figure 1 : FP amplitude field with ternary +1 / -1 symbols
%            overlaid, on a shared symmetric colour scale.
%
% Figure 2 : Gaussian-reconstructed ternary field against the FP
%            field, which is the spatial analogue of the FIR
%            reconstruction used on the time-domain probe.
%
% Reported per snapshot:
%
%   - polarity agreement, with a dead band so that cells where the
%     FP sign is meaningless are excluded rather than counted as
%     free hits
%   - the same score for randomly shuffled ternary symbols, which
%     is the chance level the real score must be judged against
%   - Pearson correlation between the reconstructed ternary field
%     and the FP field
%   - the best integer spatial shift, which should be (0,0) if the
%     two schemes propagate at the same speed
%
% Inputs: snapshot_%04d.txt and fp_snapshot_%04d.txt, as written
% by the C solver.
%
% No toolbox packages required.
% ================================================================

1;

% ------------------------------------------------
% Read one snapshot, with validation
% ------------------------------------------------
function A = read_snapshot(filename, NX, NY)

    fid = fopen(filename, 'r');

    if fid < 0
        error('Cannot open file: %s', filename);
    end

    fgetl(fid);        % "STEP n"
    fgetl(fid);        % blank

    [v, count] = fscanf(fid, '%f', [NY NX]);

    fclose(fid);

    if count ~= NX*NY
        error('%s: expected %d values, read %d.', ...
              filename, NX*NY, count);
    end

    A = v';            % A(i,j) = Ez[i][j]
end

% ------------------------------------------------
% Normalised Gaussian smoothing, edge corrected
% ------------------------------------------------
function S = gsmooth(A, sigma)

    r = max(1, ceil(3*sigma));

    [gx, gy] = meshgrid(-r:r, -r:r);

    K = exp(-(gx.^2 + gy.^2) / (2*sigma^2));
    K = K / sum(K(:));

    % Divide by the smoothed unit field so that cells near the
    % boundary are not artificially damped.
    S = conv2(A, K, 'same') ./ conv2(ones(size(A)), K, 'same');
end

% ------------------------------------------------
% Polarity agreement with a dead band
%
% Cells where |F| is below thresh*max|F| carry no usable sign and
% are excluded from both numerator and denominator.
% ------------------------------------------------
function [score, nscored, nactive] = polarity_score(T, F, thresh)

    scale = max(abs(F(:)));

    if scale == 0
        score = NaN; nscored = 0; nactive = 0;
        return;
    end

    active  = (T ~= 0);
    usable  = (abs(F) >= thresh*scale);

    scored  = active & usable;

    nactive = sum(active(:));
    nscored = sum(scored(:));

    if nscored == 0
        score = NaN;
        return;
    end

    score = 100 * sum(sign(T(scored)) == sign(F(scored))) / nscored;
end

% ------------------------------------------------
% Chance level: keep the symbol counts, destroy the
% spatial arrangement
% ------------------------------------------------
function s = chance_level(T, F, thresh, ntrial)

    sym = T(T ~= 0);
    N   = numel(T);
    acc = zeros(ntrial,1);

    for k = 1:ntrial
        R = zeros(size(T));
        idx = randperm(N, numel(sym));
        R(idx) = sym(randperm(numel(sym)));
        acc(k) = polarity_score(R, F, thresh);
    end

    s = mean(acc(~isnan(acc)));
end

% ------------------------------------------------
% Blue-white-red diverging colormap
% ------------------------------------------------
function cm = bwr(n)

    if nargin < 1, n = 256; end

    t = linspace(-1, 1, n)';

    r = min(1, max(0, 1 + t));
    b = min(1, max(0, 1 - t));
    g = 1 - abs(t);

    cm = [r g b];
end

% ------------------------------------------------
% Best integer shift by direct search
% ------------------------------------------------
function [di, dj, best] = best_shift(A, B, maxshift)

    best = -Inf; di = 0; dj = 0;

    m = maxshift + 1;
    core = @(X) X(m:end-m+1, m:end-m+1);

    b = core(B);
    b = b(:) - mean(b(:));
    nb = norm(b);

    for p = -maxshift:maxshift
        for q = -maxshift:maxshift
            a = core(circshift(A, [p q]));
            a = a(:) - mean(a(:));
            na = norm(a);
            if na == 0 || nb == 0, continue; end
            c = sum(a.*b) / (na*nb);
            if c > best
                best = c; di = p; dj = q;
            end
        end
    end
end

% ================================================================
% Main
% ================================================================

clear -x;
close all;

NX = 80;
NY = 80;

% Snapshots to show. Leave empty to use the first four found.
steps = [180 260 340 480];

FONT_SIZE      = 12;
TITLE_SIZE     = 13;

% Dead band: FP cells below this fraction of the frame maximum
% have no meaningful sign.
SIGN_THRESHOLD = 0.10;

% Width of the spatial reconstruction kernel, in cells.
SIGMA          = 2.0;

% Trials for the shuffled-symbol baseline.
NTRIAL         = 20;

% Largest spatial shift searched, in cells.
MAXSHIFT       = 4;

% ------------------------------------------------
% Discover snapshots if none were specified
% ------------------------------------------------

if isempty(steps)

    d = dir('snapshot_*.txt');
    found = [];

    for k = 1:numel(d)
        tok = regexp(d(k).name, 'snapshot_(\d+)\.txt', 'tokens');
        if ~isempty(tok)
            found(end+1) = str2double(tok{1}{1});
        end
    end

    found = sort(found);

    if numel(found) < 4
        error('Found only %d snapshots; need at least 4.', numel(found));
    end

    steps = found(round(linspace(1, numel(found), 4)));

    printf('Using snapshots: %s\n', mat2str(steps));
end

nsnap = numel(steps);

% ------------------------------------------------
% Load everything first, so the colour scale can be
% shared across panels
% ------------------------------------------------

T = cell(nsnap,1);
F = cell(nsnap,1);
R = cell(nsnap,1);
G = cell(nsnap,1);

for k = 1:nsnap

    tf = sprintf('snapshot_%04d.txt',    steps(k));
    ff = sprintf('fp_snapshot_%04d.txt', steps(k));

    T{k} = read_snapshot(tf, NX, NY);
    F{k} = read_snapshot(ff, NX, NY);

    % Spatial analogue of the time-domain FIR reconstruction.  The
    % SAME kernel is applied to the reference, exactly as the FIR is
    % applied to both traces in the time domain.  A sigma = 2 kernel
    % removes about 55 % of the signal at a 10-cell wavelength, so
    % smoothing only the ternary side would compare two different
    % band limits.
    R{k} = gsmooth(T{k}, SIGMA);
    G{k} = gsmooth(F{k}, SIGMA);
end

clim = 0;
for k = 1:nsnap
    clim = max(clim, max(abs(F{k}(:))));
end

% ------------------------------------------------
% Metrics
% ------------------------------------------------

printf('\n');
printf('=========================================================================\n');
printf(' TERNARY / FP SNAPSHOT AGREEMENT\n');
printf('=========================================================================\n');
printf('\n');
printf('  step   active   scored   polarity   chance   corr(recon,FP)   shift\n');
printf('  -----------------------------------------------------------------------\n');

score  = zeros(nsnap,1);
chance = zeros(nsnap,1);
rho    = zeros(nsnap,1);

for k = 1:nsnap

    [score(k), nscored, nactive] = ...
        polarity_score(T{k}, F{k}, SIGN_THRESHOLD);

    chance(k) = chance_level(T{k}, F{k}, SIGN_THRESHOLD, NTRIAL);

    [di, dj, rho(k)] = best_shift(R{k}, G{k}, MAXSHIFT);

    printf('  %4d   %6d   %6d   %7.1f%%  %6.1f%%   %13.4f   (%+d,%+d)\n', ...
           steps(k), nactive, nscored, score(k), chance(k), rho(k), di, dj);
end

printf('\n');
printf('  Polarity is scored only where |E_z^{fp}| >= %.2f of the frame maximum;\n', ...
       SIGN_THRESHOLD);
printf('  cells below that carry no usable sign. The chance column is the same\n');
printf('  score for randomly rearranged ternary symbols and is the level the\n');
printf('  polarity column must beat to mean anything.\n');
printf('\n');
printf('  A shift other than (0,0) means the two schemes propagate at different\n');
printf('  speeds, which is the signature of a mismatched effective CFL.\n');
printf('\n');

% ------------------------------------------------
% Figure 1: overlay
% ------------------------------------------------

f1 = figure(100);
clf;
set(f1, 'Position', [80 80 950 820]);
colormap(bwr(256));

for k = 1:nsnap

    subplot(2, 2, k);

    imagesc(F{k}, [-clim clim]);

    axis image;
    set(gca, 'YDir', 'normal', 'FontSize', FONT_SIZE);

    hold on;

    [r1, c1] = find(T{k} ==  1);
    plot(c1, r1, '.', 'MarkerSize', 5, 'Color', [0 0 0]);

    [r2, c2] = find(T{k} == -1);
    plot(c2, r2, 'o', 'MarkerSize', 2.5, 'Color', [0 0.5 0], ...
         'LineWidth', 0.5);

    hold off;

    title(sprintf('Step %d   polarity %.1f%% (chance %.1f%%)', ...
                  steps(k), score(k), chance(k)), ...
          'FontSize', TITLE_SIZE);

%    if k == 2
        h = colorbar;
        set(h, 'FontSize', FONT_SIZE);
        ylabel(h, 'E_z  (floating point)');
%    end
end

print(f1, 'snapshot_overlay.png', '-dpng', '-r120');

% ------------------------------------------------
% Figure 2: reconstruction against reference
% ------------------------------------------------

f2 = figure(101);
clf;
set(f2, 'Position', [120 120 1150 560]);
colormap(bwr(256));

for k = 1:nsnap

    % Match the reconstruction amplitude to the reference by least
    % squares, so the comparison is of shape rather than scale.
    a = R{k}(:);
    b = G{k}(:);
    g = (a'*b) / (a'*a);

    subplot(2, nsnap, k);
    imagesc(G{k}, [-clim clim]);
    axis image;
    set(gca, 'YDir', 'normal', 'FontSize', FONT_SIZE-2);
    title(sprintf('FP smoothed, step %d', steps(k)), 'FontSize', FONT_SIZE);

    subplot(2, nsnap, nsnap + k);
    imagesc(g*R{k}, [-clim clim]);
    axis image;
    set(gca, 'YDir', 'normal', 'FontSize', FONT_SIZE-2);
    title(sprintf('ternary recon, r = %.3f', rho(k)), ...
          'FontSize', FONT_SIZE);
end

print(f2, 'snapshot_reconstruction.png', '-dpng', '-r120');

printf('Written: snapshot_overlay.png, snapshot_reconstruction.png\n');
