clear; close all;

% load inputs
d = load("tmm_inputs.mat"); % wavelengths, n_p, n_m, theta_bends, n_1, n_2, theta_t1, theta_t2
[~, i1550] = min(abs(d.wavelengths - 1.55));

% define search ranges
L1 = 4:0.1:20; % straight coupling length of each symmetric coupler
L2 = 0.5:0.02:8; % straight length of the phase-control arms (between the tapers)

eta_1550 = zeros(numel(L1), numel(L2)); % cross split at 1550 nm
worst = zeros(numel(L1), numel(L2)); % largest |split - 0.5| over the C-band

for a = 1:numel(L1)
    for b = 1:numel(L2)
        eta = tmm_split(L1(a), L2(b), d);
        eta_1550(a, b) = eta(i1550);
        worst(a, b) = max(abs(eta - 0.5));
    end
end

% best: 50 +- 1 % at 1550 nm with the flattest band 
w = worst;
w(abs(eta_1550 - 0.5) > 0.01) = Inf;
[w_best, idx] = min(w(:));
[a, b] = ind2sub(size(w), idx);
fprintf("Best: L1 = %.2f um, L2 = %.2f um, split at 1550 nm = %.2f %%, worst over C-band = +-%.2f %%\n", ...
        L1(a), L2(b), 100 * eta_1550(a, b), 100 * w_best);

% ---- Maps (transposed so L1 is on the x axis and L2 on the y axis) ----
figure;
tiledlayout(1, 2);
nexttile;
imagesc(L1, L2, eta_1550.'); axis xy; colorbar; hold on;
contour(L1, L2, eta_1550.', [0.5 0.5], 'k', 'LineWidth', 2);
plot(L1(a), L2(b), 'w*', 'MarkerSize', 12);
xlabel('L_1 (\mum)'); ylabel('L_2 (\mum)'); title('Cross split at 1550 nm');
nexttile;
imagesc(L1, L2, 100 * worst.'); axis xy; colorbar; clim([0 10]); hold on;
contour(L1, L2, eta_1550.', [0.5 0.5], 'w', 'LineWidth', 2);
plot(L1(a), L2(b), 'w*', 'MarkerSize', 12);
xlabel('L_1 (\mum)'); ylabel('L_2 (\mum)'); title('Worst |split - 50| over C-band (%)');

% ---- Split vs wavelength: best TMM design and the design built in PhotonForge ----
figure; hold on;
plot(1000 * d.wavelengths, 100 * tmm_split(L1(a), L2(b), d), 'LineWidth', 2, ...
     'DisplayName', sprintf('TMM best (L_1 %.2f, L_2 %.2f)', L1(a), L2(b)));
yline([45 55], ':', 'HandleVisibility', 'off'); yline([49 51], '--', 'HandleVisibility', 'off');
xlabel('Wavelength (nm)'); ylabel('Split ratio (%)'); legend('Location', 'best'); ylim([40 60]);


function eta = tmm_split(L1, L2, d)
    % cross split ratio at every wavelength
    eta = zeros(size(d.wavelengths));
    for k = 1:numel(d.wavelengths)
        lam = d.wavelengths(k);
        % symmetric coupler: straight section plus S-bends
        th = pi * (d.n_p(k) - d.n_m(k)) / lam * L1 + d.theta_bends(k);
        C = [cos(th) 1i*sin(th); 1i*sin(th) cos(th)];
        % tapers and phase control arms
        Pt = diag(exp(1i * [d.theta_t1(k) d.theta_t2(k)]));
        P = diag(exp(1i * 2*pi / lam * [d.n_1(k) d.n_2(k)] * L2));
        E_out = C * Pt * P * Pt * C * [1; 0]; % light enters waveguide 1
        eta(k) = abs(E_out(2))^2 / (abs(E_out(1))^2 + abs(E_out(2))^2);
    end
end