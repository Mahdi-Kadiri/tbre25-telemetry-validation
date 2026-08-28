% validate_channels.m
% Sensor channel validation for AiM-logged TBRe telemetry.
% Author : Mahdi Kadiri
%
% PURPOSE
%   Establishes, from first-principles kinematics rather than channel names,
%   which logged channel is which and whether each is correctly scaled.
%   Run this before any analysis that consumes yaw rate or lateral
%   acceleration.
%
% METHOD
%   Four independent tests, each with its own reference:
%     T1  Channel identity      - correlation against V*r and dV/dt.
%                                 V*r is the kinematic lateral acceleration
%                                 and involves no accelerometer.
%     T2  Accelerometer scale   - |a| at rest must equal 1 g. Gravity is an
%                                 absolute reference, so this validates the
%                                 accelerometer triad independently of
%                                 everything else.
%     T3  Timing skew           - lag that maximises correlation between the
%                                 accelerometer and the GPS-derived
%                                 reference. Skew is per-file; never carry a
%                                 value over from another log.
%     T4  Yaw rate scale        - three routes: total least squares against
%                                 V*r, TLS against GPS course rate, and the
%                                 closed-loop heading integral. The integral
%                                 is immune to roll-gravity contamination of
%                                 the lateral channel, so it is the most
%                                 trustworthy of the three.
%
% WHY TOTAL LEAST SQUARES
%   Ordinary least squares assumes the regressor is noise-free. Both
%   quantities here are measurements, so OLS is attenuated toward zero and
%   would understate a scale error. TLS is symmetric in the two variables.
%
% LIMITATIONS
%   - Steady-state identity a_y = V*r holds only where sideslip rate is
%     small; the steady-state mask enforces this approximately, not exactly.
%   - A chassis-mounted lateral accelerometer reads g*sin(roll) in addition
%     to true lateral acceleration, inflating T4's regression routes by a
%     few percent at high lateral g. The heading integral does not suffer
%     from this.
%   - IMU offset from the CG adds x_s*rdot to measured lateral acceleration.
%     Suppressed in the steady-state set (rdot ~ 0), not elsewhere.
%   - A single static pose cannot separate sensor mounting tilt from ground
%     slope. The reported static offset is the sum of both.

clear; clc; close all;

[fn, fp] = uigetfile({'*.csv','AiM export (*.csv)'}, 'Select telemetry CSV');
if isequal(fn,0), return; end
runValidation(fullfile(fp,fn));

%% ========================================================================
function runValidation(csvPath)

G0        = 9.80665;
V_MOVING  = 5.0;      % m/s  - below this, GPS course is unreliable
V_REST    = 0.3;      % m/s  - stationary threshold for the gravity test
MAX_LAG   = 10;       % samples to search either side for timing skew
RDOT_SS   = 0.3;      % rad/s^2 - steady-state yaw acceleration limit
VDOT_SS   = 1.0;      % m/s^2   - steady-state longitudinal limit
R_MIN_SS  = 0.15;     % rad/s   - must be actually cornering

T = readtable(csvPath, 'VariableNamingRule','preserve');
[~, base] = fileparts(csvPath);
fprintf('\n===== %s =====\n', base);

t   = T.("Time");
V   = T.("speed");                    % m/s
r   = deg2rad(T.("YawRate"));         % rad/s
dt  = median(diff(t));
fprintf('%d samples at %.1f Hz, %.1f s\n', height(T), 1/dt, t(end)-t(1));

moving = V > V_MOVING;
rest   = V < V_REST;
fprintf('moving %d (%.1f%%), stationary %d\n\n', ...
        sum(moving), 100*sum(moving)/height(T), sum(rest));

% ---- GPS course and its rate (independent of the IMU) -------------------
R_EARTH = 6371000;
lat = deg2rad(T.("lat"));  lon = deg2rad(T.("lon"));
dx  = gradient(lon).*cos(lat)*R_EARTH;
dy  = gradient(lat)*R_EARTH;
course  = unwrap(atan2(dx, dy));
dcourse = gradient(course, t);

% Validate the speed channel against position differentiation before using
% it as a reference for anything else.
vGeo = hypot(dx,dy)./gradient(t);
mv   = moving & isfinite(vGeo) & vGeo < 30;
pV   = polyfit(vGeo(mv), V(mv), 1);
fprintf('[T0] speed vs lat/lon differentiation: slope %.4f, corr %.4f\n\n', ...
        pV(1), corr2v(vGeo(mv), V(mv)));

ay_kin = V .* r / G0;                 % kinematic lateral accel [g]
dV     = gradient(V, t) / G0;         % longitudinal accel [g]

%% ---- T1: channel identity ---------------------------------------------
fprintf('[T1] CHANNEL IDENTITY (no channel names trusted)\n');
cands = {'InlineAcc','LateralAcc','VerticalAcc','GXg','GYg','GZg'};
fprintf('  %-12s %8s %8s %12s %12s\n','channel','mean','std','corr V*r','corr dV/dt');
for k = 1:numel(cands)
    if ~ismember(cands{k}, T.Properties.VariableNames), continue; end
    x = T.(cands{k});
    m = moving & isfinite(x);
    fprintf('  %-12s %8.3f %8.3f %12.3f %12.3f\n', cands{k}, ...
            mean(x(m)), std(x(m)), corr2v(x(m),ay_kin(m)), corr2v(x(m),dV(m)));
end
fprintf('  -> lateral  = highest |corr| against V*r\n');
fprintf('  -> longitud = highest |corr| against dV/dt\n\n');

LAT_CH = 'InlineAcc';                 % set from T1 output
ay     = T.(LAT_CH);

%% ---- T2: accelerometer scale against gravity --------------------------
fprintf('[T2] ACCELEROMETER SCALE (gravity reference, at rest)\n');
triads = {{'InlineAcc','LateralAcc','VerticalAcc'}, {'GXg','GYg','GZg'}};
for i = 1:numel(triads)
    tri = triads{i};
    if ~all(ismember(tri, T.Properties.VariableNames)), continue; end
    A   = [T.(tri{1})(rest), T.(tri{2})(rest), T.(tri{3})(rest)];
    mag = vecnorm(A, 2, 2);
    fprintf('  %s / %s / %s\n', tri{1}, tri{2}, tri{3});
    fprintf('    static means: %+.4f %+.4f %+.4f g\n', mean(A,1));
    fprintf('    |a| at rest : %.4f g  -> scale error %+.2f%%\n', ...
            mean(mag), 100*(mean(mag)-1));
end
offset = mean(ay(rest));
fprintf('  static offset on %s: %+.4f g (= %.1f deg tilt; sensor or ground)\n\n', ...
        LAT_CH, offset, rad2deg(asin(min(abs(offset),1))));

%% ---- T3: timing skew ---------------------------------------------------
bestLag = 0; bestC = 0;
for k = -MAX_LAG:MAX_LAG
    a = shiftNaN(ay, k);
    m = moving & isfinite(a) & isfinite(ay_kin);
    c = corr2v(a(m), ay_kin(m));
    if abs(c) > abs(bestC), bestC = c; bestLag = k; end
end
m0  = moving & isfinite(ay);
fprintf('[T3] TIMING SKEW (IMU relative to GPS-derived reference)\n');
fprintf('  best lag %d samples (%.0f ms): corr %.4f  (zero-lag corr %.4f)\n', ...
        bestLag, bestLag*dt*1000, bestC, corr2v(ay(m0), ay_kin(m0)));
fprintf('  NOTE: skew is per-file. Re-derive for every log.\n\n');
ayd = shiftNaN(ay, bestLag);

%% ---- steady-state mask -------------------------------------------------
rdot = gradient(r, t);  Vdot = gradient(V, t);
ss   = moving & isfinite(ayd) & abs(rdot) < RDOT_SS & ...
       abs(Vdot) < VDOT_SS & abs(r) > R_MIN_SS;
fprintf('Steady-state samples: %d (%.1f%% of moving), corr(ay, V*r) = %.4f\n\n', ...
        sum(ss), 100*sum(ss)/sum(moving), corr2v(ayd(ss), ay_kin(ss)));

%% ---- T4: yaw rate scale ------------------------------------------------
fprintf('[T4] YAW RATE SCALE (three independent routes)\n');
g1 = tlsSlope(ayd(ss), ay_kin(ss));
fprintf('  a: TLS %s vs V*r, steady-state   -> gyro gain %.4f\n', LAT_CH, g1);

mc = moving & isfinite(dcourse) & abs(dcourse) < 3;
g2 = 1/tlsSlope(r(mc), dcourse(mc));
fprintf('  b: TLS YawRate vs GPS course rate -> gyro gain %.4f\n', g2);

% Closed-loop heading integral over contiguous moving blocks.
idx = find(moving);  brk = [0; find(diff(idx)>1); numel(idx)];
Ir = 0; Ic = 0; nb = 0;
for k = 1:numel(brk)-1
    b = idx(brk(k)+1:brk(k+1));
    if numel(b) < 200, continue; end
    Ir = Ir + trapz(t(b), r(b));
    Ic = Ic + (course(b(end)) - course(b(1)));
    nb = nb + 1;
end
g3 = Ic/Ir;
fprintf('  c: heading integral over %d blocks -> gyro gain %.4f\n', nb, g3);
fprintf('     gyro %.2f turns vs GPS %.2f turns\n', Ir/(2*pi), Ic/(2*pi));
fprintf('  -> route (c) is immune to roll-gravity error and is preferred.\n');
fprintf('  RECOMMENDED GAIN: %.3f\n\n', g3);

% Linearity: a genuine scale error is roughly constant with magnitude.
fprintf('  linearity of gyro/GPS ratio across magnitude:\n');
edges = [0.1 0.3 0.6 1.0 1.6];
for k = 1:numel(edges)-1
    s = mc & abs(dcourse) > edges(k) & abs(dcourse) <= edges(k+1);
    if sum(s) < 200, continue; end
    p = polyfit(dcourse(s), r(s), 1);
    fprintf('    |courseRate| %.1f-%.1f rad/s: n=%6d  ratio %.4f\n', ...
            edges(k), edges(k+1), sum(s), p(1));
end

%% ---- plots -------------------------------------------------------------
darkfig('Channel validation');
subplot(1,2,1); hold on; grid on;
scatter(ay_kin(m0), ay(m0), 4, [0.30 0.65 1.00], 'filled', 'MarkerFaceAlpha',0.15);
scatter(ay_kin(ss), ayd(ss), 6, [1.00 0.75 0.20], 'filled');
lim = [-2.5 2.5];
plot(lim, lim, '--', 'Color',[0.7 0.7 0.7], 'LineWidth',1.2);
plot(lim, lim*g3, '-', 'Color',[0.20 0.90 0.45], 'LineWidth',1.6);
xlabel('V\cdot r / g  [g]'); ylabel(sprintf('%s  [g]', LAT_CH));
title('Measured vs kinematic lateral acceleration');
legend({'all moving','steady-state','ideal 1:1', ...
        sprintf('gain-corrected (%.3f)', g3)}, 'Location','northwest');
styleaxes; xlim(lim); ylim(lim);

subplot(1,2,2); hold on; grid on;
scatter(dcourse(mc), rad2deg(r(mc)), 4, [0.30 0.65 1.00], 'filled', 'MarkerFaceAlpha',0.15);
xl = [-2 2];
plot(xl, rad2deg(xl), '--', 'Color',[0.7 0.7 0.7], 'LineWidth',1.2);
plot(xl, rad2deg(xl/g3), '-', 'Color',[0.20 0.90 0.45], 'LineWidth',1.6);
xlabel('GPS course rate [rad/s]'); ylabel('YawRate [deg/s]');
title(sprintf('Gyro under-reads: gain %.3f required', g3));
legend({'data','ideal 1:1','fitted gain'}, 'Location','northwest');
styleaxes; xlim(xl);

fprintf('\nApply before any downstream analysis:\n');
fprintf('  YawRate_corrected = YawRate * %.3f;\n', g3);
fprintf('  %s_corrected = %s - (%+.4f);   %% static offset\n', LAT_CH, LAT_CH, offset);
fprintf('  and deskew by %d samples (%.0f ms).\n\n', bestLag, bestLag*dt*1000);
end

%% ========================================================================
function s = tlsSlope(a, b)
% Total least squares slope of a = s*b. Symmetric in a and b, so unlike OLS
% it is not attenuated by noise in the regressor.
a = a - mean(a);  b = b - mean(b);
sxx = mean(b.^2); syy = mean(a.^2); sxy = mean(a.*b);
s = (syy - sxx + sqrt((syy-sxx)^2 + 4*sxy^2)) / (2*sxy);
end

function c = corr2v(a, b)
m = isfinite(a) & isfinite(b);
a = a(m) - mean(a(m)); b = b(m) - mean(b(m));
c = (a'*b) / (norm(a)*norm(b));
end

function y = shiftNaN(x, k)
y = nan(size(x));
if k > 0
    y(1:end-k) = x(k+1:end);
elseif k < 0
    y(-k+1:end) = x(1:end+k);
else
    y = x;
end
end

function darkfig(name)
figure('Name',name,'Color','k','Position',[80 80 1400 560]);
end

function styleaxes
ax = gca;
set(ax, 'Color','k', 'XColor','w', 'YColor','w', ...
        'GridColor',[0.45 0.45 0.45], 'GridAlpha',0.4);
set(ax.Title, 'Color','w');
lg = legend; if ~isempty(lg), set(lg,'TextColor','w','Color',[0.12 0.12 0.12]); end
end