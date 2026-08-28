% plot_yaw_vs_lateral.m
% Yaw acceleration vs lateral acceleration - measured capability envelope.
% Author : Mahdi Kadiri
% Data   : TBRe25, FSG endurance (a_3780), AiM logger, 20 Hz
%
% OBJECTIVE
%   Characterise the measured envelope of yaw acceleration against lateral
%   acceleration over a full endurance run, and quantify how yaw
%   acceleration capability varies with lateral acceleration.
%
% WHAT THIS PLOT SHOWS
%   With N = I_z * rdot, the vertical axis is yaw moment to within a
%   constant. The scatter is therefore a measured analogue of a Milliken
%   Moment Method diagram: the outer boundary is the combined lateral force
%   and yaw moment capability of the car as actually exercised on track.
%
% WHAT THIS PLOT DOES NOT SHOW
%   - It cannot attribute boundary segments to the front or rear axle. A
%     constructed MMM diagram separates axles because steer angle and body
%     slip are swept independently. Measured data gives only the convex hull
%     of whatever the driver did, with no mechanism behind any edge. An
%     earlier version of this analysis labelled edges "front tyre limit" and
%     "rear tyre limit"; those labels had no basis and have been removed.
%   - It cannot diagnose understeer or oversteer. Balance is a slip angle
%     comparison (alpha_f vs alpha_r), not a yaw acceleration observation.
%     Understeering, neutral and mildly oversteering cars all settle to
%     rdot ~ 0 in a sustained corner, so the steady-state condition is
%     agnostic to balance.
%   - The boundary is exercised, not exhaustive. It reflects driver inputs
%     over this run, so it is a lower bound on capability.
%
% SENSOR CORRECTIONS APPLIED
%   1. Yaw rate scale. The gyro under-reads. Gain established by four
%      independent routes in validate_channels.m, of which the closed-loop
%      heading integral (gyro turns vs GPS course turns over the run) is
%      preferred because it is immune to roll-gravity contamination of the
%      lateral channel.
%   2. Lateral channel identity. InlineAcc, not LateralAcc, is the true
%      lateral channel. Established by correlation against V*r (kinematic
%      lateral acceleration, no accelerometer involved), not by channel
%      name. LateralAcc is the longitudinal channel, sign inverted.
%   3. Static offset. The lateral channel reads a non-zero mean at rest,
%      equivalent to a few degrees of tilt. Sensor mounting and ground slope
%      cannot be separated from a single static pose; the sum is removed.
%   4. Timing skew between the IMU and GPS-derived quantities. Per file -
%      re-derive it for every log, never carry a value across.
%
% LIMITATIONS
%   - The IMU is offset from the CG by an unmeasured distance. Measured
%     lateral acceleration contains x_s*rdot, which is correlated with the
%     vertical axis and therefore distorts the boundary rather than adding
%     noise. Bounded by the sweep at the end of this script, not eliminated.
%   - A chassis-mounted lateral accelerometer also reads g*sin(roll), a few
%     percent at high lateral g.
%   - Yaw moment scaling depends on I_z, estimated by component summation
%     from the team mass register and uncertain to roughly +/-20%.

clear; clc; close all;

%% ---- CONFIGURATION -----------------------------------------------------
C.gyroGain   = 1.141;   % from validate_channels.m, heading-integral route
C.deskew     = -2;      % samples, IMU relative to GPS-derived quantities
C.Iz         = 100;     % kg.m^2, component summation, +/-20%
C.vMoving    = 3.0;     % m/s
C.vRest      = 0.3;     % m/s
C.smoothWin  = 5;       % samples, yaw rate smoothing before differentiation
C.binWidth   = 0.25;    % g
C.binMinN    = 50;
C.pLo        = 1;       % percentile for lower envelope
C.pHi        = 99;      % percentile for upper envelope
C.offsetSweep = [-0.4 -0.2 0 0.2 0.4];   % m, IMU longitudinal offset test

[fn, fp] = uigetfile({'*.csv','AiM export (*.csv)'}, 'Select telemetry CSV');
if isequal(fn,0), return; end
runEnvelope(fullfile(fp,fn), C);

%% ========================================================================
function runEnvelope(csvPath, C)

T = readtable(csvPath, 'VariableNamingRule','preserve');
[~, base] = fileparts(csvPath);
fprintf('\n===== %s =====\n', base);

t  = T.("Time");
V  = T.("speed");
ay = T.("InlineAcc");
rRaw = deg2rad(T.("YawRate"));

% --- corrections ---
offset = mean(ay(V < C.vRest), 'omitnan');
ay     = shiftNaN(ay - offset, C.deskew);
r      = rRaw * C.gyroGain;
fprintf('gyro gain %.3f | static offset %+.4f g | deskew %d samples (%.0f ms)\n', ...
        C.gyroGain, offset, C.deskew, C.deskew*median(diff(t))*1000);

% --- yaw acceleration ---
rdot    = gradient(movmean(r, C.smoothWin, 'omitnan'), t);
rdotRaw = gradient(movmean(rRaw, C.smoothWin, 'omitnan'), t);

m = V > C.vMoving & isfinite(ay) & isfinite(rdot);
x = ay(m);  y = rdot(m);
fprintf('samples: %d moving of %d\n\n', sum(m), height(T));

%% ---- headline numbers --------------------------------------------------
fprintf('LATERAL ACCELERATION\n');
fprintf('  peak      : %+.3f g left / %+.3f g right\n', min(x), max(x));
fprintf('  p1 / p99  : %+.3f / %+.3f g\n', pct(x,1), pct(x,99));
asym = 100*abs(abs(min(x)) - max(x)) / max(abs(min(x)), max(x));
fprintf('  asymmetry : %.1f%%  (circuit layout, not a car property -\n', asym);
fprintf('              track this across circuits rather than interpreting it)\n\n');

fprintf('YAW ACCELERATION\n');
fprintf('  range     : %+.2f .. %+.2f rad/s^2\n', min(y), max(y));
fprintf('  uncorrected would read %+.2f .. %+.2f  (%.0f%% smaller)\n', ...
        min(rdotRaw(m)), max(rdotRaw(m)), 100*(1-1/C.gyroGain));
fprintf('  peak yaw moment |N| = I_z*|rdot| = %.0f N.m at I_z = %.0f kg.m^2\n\n', ...
        C.Iz*max(abs(y)), C.Iz);

%% ---- capability taper: the substantive result --------------------------
fprintf('YAW ACCELERATION CAPABILITY vs LATERAL ACCELERATION\n');
fprintf('  %-14s %8s %14s %9s\n','|a_y| band','n','p1-p99 width','median');
bands = [0 0.25; 0.25 0.75; 0.75 1.25; 1.25 1.75; 1.75 2.5];
w = nan(size(bands,1),1);
for k = 1:size(bands,1)
    s = abs(x) >= bands(k,1) & abs(x) < bands(k,2);
    if sum(s) < C.binMinN, continue; end
    w(k) = pct(y(s),C.pHi) - pct(y(s),C.pLo);
    fprintf('  %-14s %8d %14.2f %9.2f\n', ...
        sprintf('%.2f-%.2f g',bands(k,1),bands(k,2)), sum(s), w(k), median(abs(y(s))));
end
fprintf('  -> taper %.0f%% of the moderate-g width at high lateral g.\n', 100*w(4)/w(2));
fprintf('     Mechanism: as the tyres approach lateral saturation, less\n');
fprintf('     capacity remains to generate yaw moment. This is the\n');
fprintf('     defensible reading of the envelope shape.\n\n');

%% ---- percentile envelope -----------------------------------------------
edges = -2.5:C.binWidth:2.5;
cen = []; eHi = []; eLo = [];
for k = 1:numel(edges)-1
    s = x >= edges(k) & x < edges(k+1);
    if sum(s) < C.binMinN, continue; end
    cen(end+1) = mean(edges(k:k+1)); %#ok<AGROW>
    eHi(end+1) = pct(y(s), C.pHi); %#ok<AGROW>
    eLo(end+1) = pct(y(s), C.pLo); %#ok<AGROW>
end

%% ---- IMU offset sensitivity -------------------------------------------
% Measured lateral acceleration at the IMU is a_y,CG + x_s*rdot. x_s is
% unmeasured, so bound its effect rather than ignoring it.
fprintf('IMU LONGITUDINAL OFFSET SENSITIVITY (x_s unmeasured)\n');
fprintf('  %8s %14s %14s\n','x_s [m]','peak |a_y| [g]','env width [g]');
for xs = C.offsetSweep
    xc = x - xs*y/9.80665;
    fprintf('  %8.2f %14.3f %14.3f\n', xs, max(abs(xc)), pct(xc,99)-pct(xc,1));
end
fprintf('  -> report this band until the mounting position is measured.\n\n');

%% ---- plot --------------------------------------------------------------
figure('Color','k','Position',[80 80 1000 640]); hold on; grid on;
scatter(x, y, 5, [0.30 0.65 1.00], 'filled', 'MarkerFaceAlpha', 0.18);
plot(cen, eHi, '-o', 'Color',[1.00 0.75 0.20], 'LineWidth',2, ...
     'MarkerFaceColor',[1.00 0.75 0.20], 'MarkerSize',4);
plot(cen, eLo, '-o', 'Color',[1.00 0.75 0.20], 'LineWidth',2, ...
     'MarkerFaceColor',[1.00 0.75 0.20], 'MarkerSize',4);
xlabel('Lateral acceleration [g]   (+ve = right)');
ylabel('Yaw acceleration [rad/s^2]');
title({sprintf('Measured yaw acceleration envelope - %s', strrep(base,'_',' ')), ...
       sprintf('gyro gain %.3f applied | envelope = %d-%d percentile in %.2f g bins', ...
               C.gyroGain, C.pLo, C.pHi, C.binWidth)});
legend({'moving samples', sprintf('%d/%d percentile envelope',C.pLo,C.pHi)}, ...
       'Location','northeast');

% Capture the LEFT axis limits BEFORE switching sides: after yyaxis right,
% ylim() returns the right axis's own (default) limits, not the left's.
ylLeft = ylim;
yyaxis right
ylim(ylLeft*C.Iz);
ylabel(sprintf('Yaw moment [N.m]  (I_z = %.0f kg.m^2, \\pm20%%)', C.Iz));

ax = gca;
set(ax,'Color','k','XColor','w','GridColor',[0.45 0.45 0.45],'GridAlpha',0.4);
ax.YAxis(1).Color = 'w';  ax.YAxis(2).Color = [0.6 0.6 0.6];
set(ax.Title,'Color','w');
lg = legend; set(lg,'TextColor','w','Color',[0.12 0.12 0.12]);
end

%% ========================================================================
function p = pct(x, q)
% Percentile without the Statistics Toolbox (linear interpolation, matches
% prctile for the sizes used here).
x = sort(x(isfinite(x)));
n = numel(x);
if n == 0, p = NaN; return; end
idx = 1 + (q/100)*(n-1);
lo = floor(idx); hi = ceil(idx);
p  = x(lo) + (idx-lo).*(x(hi)-x(lo));
end

function y = shiftNaN(x, k)
y = nan(size(x));
if k > 0,     y(1:end-k) = x(k+1:end);
elseif k < 0, y(-k+1:end) = x(1:end+k);
else,         y = x;
end
end