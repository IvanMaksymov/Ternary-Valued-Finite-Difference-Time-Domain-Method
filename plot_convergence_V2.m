% ================================================================
% plot_convergence.m
%
% GNU Octave script that plots the convergence study and prints the
% manuscript table.
%
% Panel 1 : relative field error against CFL, log-log, with an
%           O(CFL) reference line.
% Panel 2 : ratio of ternary to floating-point field energy, which
%           converges to unity from above.
% Panel 3 : a cut through row i = 40 of the Ez field at equal
%           physical time, for three values of CFL, against the
%           floating-point reference.
%
% The convergence numbers are READ FROM FILE, not hardcoded, so
% that the figure and the manuscript table cannot drift apart.
%
% Required inputs, all written by ternary_convergence_test.c:
%
%   convergence.txt   CFL  steps  NRMSE  E_tern/E_fp  saturations
%   cut_10000.txt     CFL = 0.1        (j, ternary_Ez, fp_Ez)
%   cut_2500.txt      CFL = 0.025
%   cut_625.txt       CFL = 0.00625
%
% Generate them with:
%
%   rm -f convergence.txt
%   for p in "10000 2000" "5000 4000" "2500 8000" \
%            "1250 16000" "625 32000"; do
%     set -- $p
%     gcc -O2 -DMC=$1 -DNSTEP=$2 -o cv ternary_convergence_test.c -lm
%     ./cv
%   done
%
% convergence.txt is appended to, so a repeated sweep adds rows;
% this script keeps the last entry for each distinct CFL.
%
% Usage:  octave plot_convergence.m
%
% No toolbox packages are required.
% ================================================================

clear;
close all;

DATAFILE = 'convergence.txt';

% ------------------------------------------------
% Read the convergence data
% ------------------------------------------------

if exist(DATAFILE, 'file') == 2

    raw = dlmread(DATAFILE);

    if columns(raw) < 4
        error('%s must have at least four columns.', DATAFILE);
    end

    % Keep the last row for each distinct CFL, then sort descending
    % so the table reads coarse to fine.
    [u, ~, ic] = unique(raw(:,1));

    keep = zeros(numel(u), 1);
    for k = 1:numel(u)
        rows = find(ic == k);
        keep(k) = rows(end);
    end

    raw = raw(keep, :);
    [~, order] = sort(raw(:,1), 'descend');
    raw = raw(order, :);

else

    warning('%s not found; falling back to built-in values.', DATAFILE);
    
    STOP

%    raw = [0.10000  2000  1.943712  4.309842  0;
%           0.05000  4000  1.013461  1.658087  0;
%           0.02500  8000  0.516192  1.101137  0;
%           0.01250 16000  0.280440  1.001824  0;
%           0.00625 32000  0.145070  0.983370  0];

end

cfl    = raw(:,1);
steps  = raw(:,2);
err    = raw(:,3);
eratio = raw(:,4);

if numel(cfl) < 2
    error('Need at least two runs to plot a convergence study.');
end

% ------------------------------------------------
% Fitted order and manuscript table
% ------------------------------------------------

p = polyfit(log(cfl), log(err), 1);

printf('\n');
printf('=================================================================\n');
printf(' CONVERGENCE AT FIXED PHYSICAL TIME\n');
printf('=================================================================\n');
printf('\n');
printf('   CFL        Steps   Relative error   E_tern/E_fp\n');
printf('   ---------------------------------------------------\n');

for k = 1:numel(cfl)
    printf('   %8.5f  %7d   %12.4f   %11.4f\n', ...
           cfl(k), steps(k), err(k), eratio(k));
end

printf('\n');
printf('   Fitted convergence order : %.4f\n', p(1));
printf('   Energy ratio at finest   : %.4f\n', eratio(end));
printf('\n');

% ------------------------------------------------
% The same rows, ready to paste into the manuscript
% ------------------------------------------------

printf('LaTeX table body:\n\n');

for k = 1:numel(cfl)
    printf('$%g$ & $%d$ & $%.4f$ & $%.4f$ \\\\\\\\\n', ...
           cfl(k), steps(k), err(k), eratio(k));
end

printf('\n');

% ------------------------------------------------
% Field cuts
% ------------------------------------------------

cut_files  = {'cut_10000.txt', 'cut_2500.txt', 'cut_625.txt'};
cut_labels = {'ternary, CFL=0.1', ...
              'ternary, CFL=0.025', ...
              'ternary, CFL=0.00625'};
cut_colors = [0.90 0.10 0.10;   % red
              0.00 0.60 0.20;   % green
              0.00 0.35 0.90];  % blue

have_cuts = true;
for k = 1:numel(cut_files)
    if exist(cut_files{k}, 'file') ~= 2
        warning('Missing %s; skipping the field-cut panel.', cut_files{k});
        have_cuts = false;
    end
end

% ------------------------------------------------
% Figure
% ------------------------------------------------

if have_cuts
    npanel = 3;
    fig = figure('Position', [100 100 1400 420]);
else
    npanel = 2;
    fig = figure('Position', [100 100 950 420]);
end

% ---------- Panel 1: error ----------

subplot(1, npanel, 1);

loglog(cfl, err, 'o-', ...
       'LineWidth', 1.6, 'MarkerSize', 6, ...
       'MarkerFaceColor', [0.12 0.47 0.71], ...
       'Color', [0.12 0.47 0.71]);
hold on;

% Reference line anchored to the coarsest point, so it is a genuine
% comparison rather than a fitted curve.
ref = err(1) * (cfl / cfl(1));
loglog(cfl, ref, 'k--', 'LineWidth', 1.0);
hold off;

grid on;
%xlabel('CFL (= 1/oversampling ratio)');
xlabel('Courant number S');
ylabel('relative field error at fixed physical time');
title(sprintf('Convergence, fitted order %.2f', p(1)));
legend('ternary vs FP  (NRMSE)', 'O(CFL) reference', ...
       'Location', 'northwest');
axis tight;

% ---------- Panel 2: energy ratio ----------

subplot(1, npanel, 2);

semilogx(cfl, eratio, 's-', ...
         'LineWidth', 1.6, 'MarkerSize', 6, ...
         'MarkerFaceColor', [0.47 0.20 0.55], ...
         'Color', [0.47 0.20 0.55]);
hold on;

xl = [min(cfl)/1.5, max(cfl)*1.5];
plot(xl, [1 1], 'k--', 'LineWidth', 1.0);
xlim(xl);
hold off;

grid on;
%xlabel('CFL (= 1/oversampling ratio)');
xlabel('Courant number S');
ylabel('E_{ternary} / E_{floating point}');
title('Field energy ratio');
legend('measured', 'exact equivalence', 'Location', 'northeast');

% ---------- Panel 3: field cuts ----------

if have_cuts

    subplot(1, npanel, 3);
    hold on;

    for k = 1:numel(cut_files)
        d = dlmread(cut_files{k});
        plot(d(:,1), d(:,2), 'LineWidth', 1.0, 'Color', cut_colors(k,:));
    end

    d = dlmread(cut_files{end});
    plot(d(:,1), d(:,3), 'k:', 'LineWidth', 2.0);
    hold off;

    grid on;
    xlabel('j  (row i=40)');
    ylabel('E_z');
    title('Field cut at equal physical time');
    legend(cut_labels{1}, cut_labels{2}, cut_labels{3}, ...
           'floating-point reference', 'Location', 'northeast');
    box on;

end

% ------------------------------------------------
% Save
% ------------------------------------------------

print(fig, 'convergence.png', '-dpng', '-r130');
%print(fig, 'convergence.pdf', '-dpdf', '-painters');
printf('Written: convergence.png\n');
