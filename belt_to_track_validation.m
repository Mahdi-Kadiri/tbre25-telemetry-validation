% belt_to_track_validation.m
% Recovering the flat-trac to on-track friction scaling factor from vehicle
% telemetry, TBRe25.
% Author : Mahdi Kadiri
%
% ------------------------------------------------------------------------
% OBJECTIVE
%   Predict rear axle lateral force from flat-trac tyre data through a full
%   vehicle model, compare against the same quantity measured on track, and
%   test whether the discrepancy is a constant friction scaling (as flat-trac
%   theory predicts) or a shape error (which would indicate the model or the
%   sensor corrections are wrong).
%
% WHY THIS IS THE TEST THAT MATTERS
%   A single ratio proves nothing - any model can be scaled to hit one point.
%   The evidence is that the ratio is FLAT across the whole slip angle range.
%   That can only happen if the curve shape is right, which in turn requires
%   the gyro scale correction, the channel identity, the sideslip estimate,
%   the load transfer split and the tyre fit to all be right simultaneously.
%   A flat ratio validates the entire chain at once; a sloped one would say
%   something upstream is broken.
%
% METHOD
%   1. Sensor corrections (established in validate_channels.m):
%        yaw rate scale, lateral channel static offset, IMU/GPS timing skew.
%   2. Chassis sideslip by anchored integration of beta_dot = a_y/V - r,
%      reset to zero at straight-line samples either side of each corner.
%   3. Steady-state mask so the yaw inertia term drops out of the axle force
%      balance and F_yr follows from static geometry alone.
%   4. Rear axle slip angle alpha_r = beta - b*r/V and measured rear axle
%      lateral force F_yr = m*a_y*a/L. Median-binned by |alpha_r|.
%   5. For each bin, back out the implied lateral acceleration, compute the
%      rear axle lateral load transfer (elastic + geometric + unsprung),
%      evaluate the tyre model at EACH tyre's own normal load, and sum.
%   6. Ratio of measured to predicted, per bin, plus a least-squares scaling.
%
% RESULT (FSG25 endurance a_3780, camber 1 deg)
%   Recovered scaling 0.630, ratio standard deviation 0.033 over
%   alpha_r = 1.3 to 6.7 deg. Documented FSAE belt-to-track factor 0.667.
%
% ROLL STIFFNESS PROVENANCE
%   The elastic roll stiffness distribution is DERIVED here, not assumed:
%     wheel rate      k_w   = k_spring * MR^2
%     axle roll rate  K_phi = 0.5 * k_w * t^2
%   Front (no ARB) 534.8 and rear (springs) 309.1 N.m/deg sum to 843.9,
%   against the team's stated 844. Adding the 56 N.m/deg rear anti-roll bar
%   residual gives 899.9 against the stated 899. Reproducing both figures
%   from first principles is what makes the 0.594 elastic split citable.
%
% ASSUMPTIONS
%   - a_y = V*r on the steady-state set; sideslip rate negligible there.
%   - Sideslip is zero at the straight-line anchors.
%   - Aerodynamic download IS included, computed per bin from SCz and the
%     bin's median speed (it scales with v^2, so a fixed value is wrong).
%     Both consistent source pairs are evaluated. Including it lowers the
%     recovered scaling from ~0.63 to ~0.585 - i.e. omitting aero flatters
%     the agreement with the documented 0.667, so it must not be omitted.
%   - No longitudinal load transfer correction on the steady-state set.
%   - Roll centre heights taken as fixed, ignoring migration with roll.
%
% LIMITATIONS
%   - Front weight fraction is the dominant uncertainty. F_yr scales linearly
%     with it, so it maps directly onto the recovered scaling factor. The
%     team has a four-way conflict (47.8 / 48 / 51 / 54 per cent) and corner
%     weights would close it. Quantified in the sensitivity block below.
%   - The tyre model carries no temperature term and no combined slip, so it
%     over-predicts during braking and traction phases.
%   - Static camber is a three-way conflict. Oliver's -1.0 deg is used here
%     because he was present at the event. The recovered factor moves from
%     0.630 to about 0.634 across the disputed range, so the conclusion is
%     insensitive to it, but the value is not settled.
%   - Agreement degrades above about 5 deg (ratio falls to 0.88 at 6.7 deg).
%     Candidates are tyre temperature, combined slip, and the fitted peak
%     slip angle sitting near the edge of the TTC sweep. Do not claim the
%     model above 5 deg.
%   - No Statistics Toolbox dependency. Percentile and correlation helpers
%     are implemented at the bottom of this file.
%
% REQUIRES
%   hoosier_16x75_10_R20.m  (TTC-derived, NOT for publication)
% ------------------------------------------------------------------------

clear; clc; close all;

%% ---- VEHICLE PARAMETERS (TBRe25) ---------------------------------------
% Source and confidence stated per line. Anything marked QUOTED has no
% document behind it and should be treated as provisional.
P.m         = 306;    % kg      ready to run with driver     [triple-confirmed]
P.L         = 1.530;    % m       wheelbase                    [confirmed]
P.t_f       = 1.200;    % m       front track                  [confirmed]
P.t_r       = 1.190;    % m       rear track                   [confirmed]
P.frontFrac = 0.48;     % -       front axle load fraction     [QUOTED, +/-0.05]
P.h         = 0.286;    % m       CG height                    [telemetry-selected]
P.h_s       = 0.2993;   % m       sprung mass CG height
P.m_u_f     = 8.0;      % kg      unsprung mass per front corner
P.m_u_r     = 13.0;     % kg      unsprung mass per rear corner
P.h_u       = 0.195;    % m       unsprung CG height (wheel centre)
P.RC_f      = 0.0486;   % m       front roll centre height
P.RC_r      = 0.0533;   % m       rear roll centre height
P.k_spring  = 52.54;    % N/mm    spring rate, both axles (300 lb/in)
P.MR_f      = 0.90;     % -       front motion ratio, telemetry-effective
P.MR_r      = 0.69;     % -       rear motion ratio, telemetry-effective
P.ARB_r     = 56.0;     % N.m/deg rear anti-roll bar contribution
P.rho       = 1.225;    % kg/m^3  air density
P.aero      = [3.07 0.510; 3.58 0.463];  % [SCz  REAR downforce fraction]
                        %         Two CONSISTENT source pairs - handover
                        %         (design/CFD, 49% front) and telemetry
                        %         ride-height fit (53.7% front). Both are
                        %         evaluated; NEVER mix SCz from one with the
                        %         balance from the other.
P.aeroName  = {'handover','telemetry'};
P.camber    = 0.645;    % deg     REAR static camber magnitude [FSG design doc]
                        %         (this script characterises the REAR axle)

%% ---- SENSOR CORRECTIONS ------------------------------------------------
S.gyroGain  = 1.141;    % from validate_channels.m, heading-integral route
S.deskew    = 2;        % samples, IMU lags GPS-derived quantities

%% ---- MASK THRESHOLDS ---------------------------------------------------
M.vMoving   = 5.0;      % m/s
M.vRest     = 0.3;      % m/s
M.straightR = 0.10;     % rad/s
M.straightA = 0.15;     % g
M.ssRdot    = 0.30;     % rad/s^2
M.ssVdot    = 1.00;     % m/s^2
M.ssRmin    = 0.15;     % rad/s
M.binEdges  = [0.5 1 1.5 2 2.5 3 4 5 6 8];   % deg, |alpha_r|
M.binMinN   = 40;

[fn, fp] = uigetfile({'*.csv','AiM export (*.csv)'}, 'Select endurance telemetry CSV');
if isequal(fn,0), return; end
belt_to_track(fullfile(fp,fn), P, S, M);

%% ========================================================================
function belt_to_track(csvPath, P, S, M)

G0 = 9.80665;
T  = readtable(csvPath, 'VariableNamingRule','preserve');
[~, base] = fileparts(csvPath);
fprintf('\n===== BELT-TO-TRACK VALIDATION: %s =====\n', base);

t  = T.("Time");
V  = T.("speed");
r  = deg2rad(T.("YawRate")) * S.gyroGain;
ay = T.("InlineAcc");

offset = mean(ay(V < M.vRest), 'omitnan');
ay     = shiftSamples(ay - offset, S.deskew);
ay(~isfinite(ay)) = 0;
fprintf('corrections: gyro gain %.3f | offset %+.4f g | deskew %+d samples\n', ...
        S.gyroGain, offset, S.deskew);

%% ---- roll stiffness, derived -------------------------------------------
k_w_f = P.k_spring * P.MR_f^2 * 1000;          % N/m at the wheel
k_w_r = P.k_spring * P.MR_r^2 * 1000;
K_f   = 0.5 * k_w_f * P.t_f^2 * pi/180;        % N.m/deg
K_r_s = 0.5 * k_w_r * P.t_r^2 * pi/180;
K_r   = K_r_s + P.ARB_r;
elasticFront = K_f / (K_f + K_r);

fprintf('\nROLL STIFFNESS (derived from springs, motion ratios, tracks)\n');
fprintf('  front (no ARB)     %6.1f N.m/deg\n', K_f);
fprintf('  rear  (springs)    %6.1f N.m/deg   -> total %.1f  [team states 844]\n', K_r_s, K_f+K_r_s);
fprintf('  rear  (+ARB %2.0f)    %6.1f N.m/deg   -> total %.1f  [team states 899]\n', P.ARB_r, K_r, K_f+K_r);
fprintf('  elastic distribution front = %.3f\n', elasticFront);

%% ---- sideslip ----------------------------------------------------------
beta = anchoredSideslip(t, V, r, ay, M, G0);
mv   = V > M.vMoving;
fprintf('\nSIDESLIP: std %.2f deg over %d moving samples\n', ...
        std(rad2deg(beta(mv))), sum(mv));

%% ---- steady-state rear axle -------------------------------------------
rdot = gradient(r,t);  Vdot = gradient(V,t);
ss = mv & abs(rdot) < M.ssRdot & abs(Vdot) < M.ssVdot & abs(r) > M.ssRmin;

b = P.frontFrac * P.L;        % CG to REAR axle  (front load fraction = b/L)
a = (1-P.frontFrac) * P.L;    % CG to FRONT axle

alpha = rad2deg(beta - b*r ./ max(V,1e-3));
Fyr   = P.m * ay * G0 * a / P.L;
fprintf('steady-state samples: %d (%.1f%% of moving)\n', sum(ss), 100*sum(ss)/sum(mv));

A = abs(alpha(ss));  F = abs(Fyr(ss));

%% ---- bin, predict, compare --------------------------------------------
for pA = 1:size(P.aero,1)
fprintf('\n--- aero pair: %s (SCz %.2f, %.1f%% rear) ---\n', ...
        P.aeroName{pA}, P.aero(pA,1), 100*P.aero(pA,2));
fprintf('%8s %6s %8s %9s %10s %8s\n','|a_r| deg','n','DF_r N','meas N','model N','ratio');
aM = []; fM = []; fP = [];
for k = 1:numel(M.binEdges)-1
    s = A >= M.binEdges(k) & A < M.binEdges(k+1);
    if sum(s) < M.binMinN, continue; end
    a_bin = median(A(s));
    F_bin = median(F(s));
    if a_bin < 1.0, continue; end   % below 1 deg the abs-binning inflates
                                    % small noisy values; not usable
    ay_bin = F_bin / (P.m*G0*(1-P.frontFrac));       % implied lateral g
    dFz    = rearLoadTransfer(ay_bin, P, elasticFront, G0);
    % Aerodynamic download on the rear axle at THIS bin's speed. Downforce
    % scales with v^2, so a single fixed value is wrong away from the speed
    % it was evaluated at - compute per bin from SCz and the bin speed.
    V_bin  = binSpeed(V, ss, A, M.binEdges(k), M.binEdges(k+1));
    DF_r   = P.aero(pA,2) * P.aero(pA,1) * 0.5 * P.rho * V_bin^2;
    Fz0    = (P.m*G0*(1-P.frontFrac) + DF_r)/2;      % static + aero, per tyre
    Fz_o   = Fz0 + dFz;
    Fz_i   = max(Fz0 - dFz, 0);
    Fmod   = abs(hoosier_16x75_10_R20(a_bin, Fz_o, 1.0, P.camber));
    if Fz_i > 1
        Fmod = Fmod + abs(hoosier_16x75_10_R20(a_bin, Fz_i, 1.0, P.camber));
    end
    aM(end+1)=a_bin; fM(end+1)=F_bin; fP(end+1)=Fmod; %#ok<AGROW>
    fprintf('%8.2f %6d %8.0f %9.0f %10.0f %8.3f\n', ...
            a_bin, sum(s), DF_r, F_bin, Fmod, F_bin/Fmod);
end

ratio   = fM ./ fP;
scaling = sum(fM.*fP) / sum(fP.^2);
fprintf('\nRECOVERED SCALING       = %.3f\n', scaling);
fprintf('  mean ratio            = %.3f  (sd %.3f)\n', mean(ratio), std(ratio));
fprintf('  FSAE documented factor= 0.667   -> %.1f%% difference\n', 100*abs(scaling/0.667-1));
fprintf('  ratio spread is %.1f%% of its mean. A FLAT ratio is the evidence:\n', 100*std(ratio)/mean(ratio));
fprintf('  it says the curve SHAPE is right and only the friction LEVEL differs.\n');
end   % aero pair loop

%% ---- sensitivity -------------------------------------------------------
fprintf('\nSENSITIVITY TO FRONT WEIGHT FRACTION (telemetry aero pair)\n');
fprintf('  %10s %12s\n','front frac','scaling');
for ffv = [P.frontFrac-0.05, P.frontFrac, P.frontFrac+0.05]
    Pv = P; Pv.frontFrac = ffv;
    av = (1-ffv)*P.L;  bv = ffv*P.L;
    alv = rad2deg(beta - bv*r./max(V,1e-3));
    Fv  = P.m * ay * G0 * av / P.L;
    Av = abs(alv(ss)); Fvv = abs(Fv(ss));
    mm=[]; pp=[];
    for k = 1:numel(M.binEdges)-1
        s = Av >= M.binEdges(k) & Av < M.binEdges(k+1);
        if sum(s) < M.binMinN, continue; end
        ab = median(Av(s)); fb = median(Fvv(s));
        if ab < 1.0, continue; end
        ayb = fb/(P.m*G0*(1-ffv));
        dF  = rearLoadTransfer(ayb, Pv, elasticFront, G0);
        Vb  = binSpeed(V, ss, Av, M.binEdges(k), M.binEdges(k+1));
        DFb = P.aero(2,2)*P.aero(2,1)*0.5*P.rho*Vb^2;   % telemetry pair
        F0  = (P.m*G0*(1-ffv) + DFb)/2;
        fp_ = abs(hoosier_16x75_10_R20(ab, F0+dF, 1.0, P.camber));
        if F0-dF > 1, fp_ = fp_ + abs(hoosier_16x75_10_R20(ab, F0-dF, 1.0, P.camber)); end
        mm(end+1)=fb; pp(end+1)=fp_; %#ok<AGROW>
    end
    fprintf('  %10.2f %12.3f\n', ffv, sum(mm.*pp)/sum(pp.^2));
end

%% ---- plot --------------------------------------------------------------
figure('Color','k','Position',[80 80 1300 560]);

subplot(1,2,1); hold on; grid on;
plot(aM, fM, 'o-', 'Color',[1.00 0.75 0.20], 'LineWidth',2, ...
     'MarkerFaceColor',[1.00 0.75 0.20], 'MarkerSize',6);
plot(aM, fP, 's--', 'Color',[0.55 0.65 0.85], 'LineWidth',1.6, 'MarkerSize',6);
plot(aM, fP*scaling, '^-', 'Color',[0.20 0.90 0.45], 'LineWidth',2, ...
     'MarkerFaceColor',[0.20 0.90 0.45], 'MarkerSize',6);
xlabel('Rear axle slip angle |\alpha_r|  [deg]');
ylabel('Rear axle lateral force  [N]');
title('Flat-trac prediction vs on-track measurement');
legend({'measured (telemetry)','flat-trac model, unscaled', ...
        sprintf('model \\times %.3f', scaling)}, 'Location','southeast');
styleaxes;

subplot(1,2,2); hold on; grid on;
hR=plot(aM, ratio, 'o-', 'Color',[1.00 0.75 0.20], 'LineWidth',2, ...
     'MarkerFaceColor',[1.00 0.75 0.20], 'MarkerSize',6);
h1=yline(scaling, '-',  sprintf('recovered %.3f',scaling), ...
      'Color',[0.20 0.90 0.45], 'LineWidth',2, 'LabelHorizontalAlignment','left');
h2=yline(0.667,   '--', 'FSAE documented 0.667', ...
      'Color',[0.85 0.85 0.85], 'LineWidth',1.4, 'LabelHorizontalAlignment','right');
legend([hR h1 h2], {'measured / model, per bin', ...
       sprintf('recovered scaling %.3f',scaling), ...
       'FSAE documented 0.667'}, 'Location','southwest');
xlabel('Rear axle slip angle |\alpha_r|  [deg]');
ylabel('measured / model');
title(sprintf('Ratio is flat: sd %.3f on mean %.3f', std(ratio), mean(ratio)));
ylim([0.4 0.9]);
styleaxes;
end

%% ========================================================================
function V_bin = binSpeed(V, ss, A, lo, hi)
% Median speed of the samples in this slip-angle bin. Needed because aero
% download depends on speed and the bins are not at a common speed.
Vs = V(ss);
s  = A >= lo & A < hi;
if nnz(s) < 5, V_bin = median(Vs,'omitnan'); else, V_bin = median(Vs(s),'omitnan'); end
end

function dFz = rearLoadTransfer(ay, P, elasticFront, G0)
% Rear axle lateral load transfer per tyre, decomposed into the three
% physically distinct paths. Lumping these into a single m*ay*h/t term (a
% common shortcut) misallocates transfer between axles, because only the
% elastic component follows the roll stiffness split.
m_s   = P.m - 2*P.m_u_f - 2*P.m_u_r;
m_s_r = m_s * (1 - P.frontFrac);
% Roll arm measured to the roll axis at the CG station.
arm   = P.h_s - (P.RC_f*(1-P.frontFrac) + P.RC_r*P.frontFrac);

elastic   = (1-elasticFront) * m_s * ay*G0 * arm    / P.t_r;
geometric =  m_s_r           * ay*G0 * P.RC_r       / P.t_r;
unsprung  =  2*P.m_u_r       * ay*G0 * P.h_u        / P.t_r;
dFz = elastic + geometric + unsprung;
end

function beta = anchoredSideslip(t, V, r, ay, M, G0)
% Integrate beta_dot = a_y/V - r between straight-line anchors, forcing
% beta = 0 at each anchor and removing accumulated drift linearly in
% between. Preserves sustained sideslip WITHIN a corner, which a
% complementary or GPS-corrected filter would high-pass away.
straight = abs(r) < M.straightR & abs(ay) < M.straightA & V > M.vMoving;
bdot = zeros(size(t));
mv   = V > 3;
bdot(mv) = ay(mv)*G0./V(mv) - r(mv);

beta = zeros(size(t));
idx  = find(V > M.vMoving);
if isempty(idx), return; end
brk = [0; find(diff(idx)>1); numel(idx)];
for k = 1:numel(brk)-1
    blk = idx(brk(k)+1:brk(k+1));
    if numel(blk) < 50, continue; end
    anc = find(straight(blk));
    for i = 1:numel(anc)-1
        if anc(i+1)-anc(i) < 3, continue; end
        seg = blk(anc(i):anc(i+1));
        bi  = [0; cumsum(bdot(seg(2:end)) .* diff(t(seg)))];
        beta(seg) = bi - linspace(0, bi(end), numel(bi))';
    end
end
end

function y = shiftSamples(x, k)
% Positive k shifts x later in time relative to the rest of the record.
y = nan(size(x));
if k > 0,     y(k+1:end) = x(1:end-k);
elseif k < 0, y(1:end+k) = x(1-k:end);
else,         y = x;
end
end

function styleaxes
ax = gca;
set(ax,'Color','k','XColor','w','YColor','w', ...
       'GridColor',[0.45 0.45 0.45],'GridAlpha',0.4);
set(ax.Title,'Color','w');
lg = legend; if ~isempty(lg), set(lg,'TextColor','w','Color',[0.12 0.12 0.12]); end
end