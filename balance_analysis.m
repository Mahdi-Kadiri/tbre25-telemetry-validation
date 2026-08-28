% balance_analysis.m
% TBRe25 steady-state balance and limiting-axle analysis by tyre-model
% inversion - no steering-angle channel required.
% Author : Mahdi Kadiri
%
% ------------------------------------------------------------------------
% OBJECTIVE
%   The car logs no steering angle, so front axle slip angle - and therefore
%   the classical understeer gradient - is not directly measurable. This
%   recovers it by inversion, and asks which axle limits peak lateral
%   acceleration, with every uncertain input bracketed rather than chosen.
%
% METHOD
%   In steady state both axle lateral forces follow from geometry alone:
%       F_yf = m*ay*g*b/L = m*ay*g*ff       F_yr = m*ay*g*(1-ff)
%   Given a TTC-derived tyre model whose scaling was recovered against
%   measured rear-axle force (belt_to_track_validation.m), INVERT it: find
%   the slip angle each axle needs to deliver its requirement at its own
%   per-tyre loads, including lateral load transfer, speed-dependent aero
%   download, per-axle camber and front toe-out.
%
%   Two metrics are reported:
%     (a) alpha_f - alpha_r  -> sub-limit balance (understeer if positive)
%     (b) AXLE UTILISATION   -> required force / axle capacity, per axle.
%         The axle reaching 100% first is the limiting axle. The DIFFERENCE
%         between the two is the margin, and is the honest way to state a
%         limiting-axle result: a 1% difference is not a verdict.
%
%   VALIDATION STATUS - the two halves are NOT equally trustworthy:
%     REAR  : validated. Model alpha_r matches telemetry-derived alpha_r to
%             ~0.2 deg over 0.5-1.25 g.
%     FRONT : PREDICTION. No measurable alpha_f exists on this car. Same
%             tyre, but extrapolated to different load, camber and toe.
%
% HEADLINE FINDING
%   Across the full uncertainty space the two axles sit within 0.2-4% of
%   each other in utilisation at the limit. The car is close to
%   SIMULTANEOUSLY LIMITED, and the limiting axle is INDETERMINATE on this
%   parameter set - it flips with both the front weight fraction and the
%   choice of aero source. That is the correct conclusion, not a failure to
%   reach one: a 1% margin difference cannot support an axle verdict.
%
%   Sub-limit balance is consistently but weakly oversteer-biased
%   (alpha_f - alpha_r = -0.19 to -0.25 deg at 1.0 g), driven mainly by the
%   1 deg front toe-out.
%
%   Predicted limit 1.43-1.46 g. Measured sustained maxima are a LOWER
%   BOUND on capability, so only bands where measurement EXCEEDS the model
%   are binding: the car demonstrably reached 1.589 g where the model says
%   1.464, so the model under-predicts by at least 9%. Bands below the
%   model are uninformative, not agreement. Matching the measured peak
%   needs
%   scaling ~0.62-0.64 against the 0.585 fitted at mid-range slip. Likely
%   cause: endurance mid-range data carries combined-slip contamination
%   (braking and traction consuming lateral capacity), suppressing measured
%   force and dragging the mid-range fit low. Plausible, unproven - stated
%   as an open discrepancy, not explained away.
%
% AERO - two CONSISTENT source pairs. NEVER MIX THEM.
%   handover   SCz 3.07, 49.0% front   (design/CFD, fixed 40 mm, 18 m/s)
%   telemetry  SCz 3.58, 53.7% front   (ride-height fit, real RH + squash)
%   Rear download is similar between pairs, so the recovered scaling is
%   pair-insensitive (0.586 vs 0.584). Front download is not, which is why
%   the limiting axle flips.
%
% OTHER STATED UNCERTAINTIES
%   front weight fraction : 0.48 quoted, four-way document conflict
%                           (47.8/48/51/54) -> corner weights would close it
%   toe convention        : datasheet says "1 deg toe out front", total vs
%                           per-side unstated. FS sheets conventionally
%                           quote TOTAL, and 1 deg/side would be an extreme
%                           drag penalty, so 0.5 deg/side is taken as
%                           likely. Affects magnitude, not sign, and drops
%                           out of the limit entirely (toe changes the slip
%                           the front needs, not its peak capacity).
%   roll centres          : 48.6/53.3 used; 23/34 and 17/30 also exist in
%                           team documents. Affects the geometric transfer
%                           term only.
%   camber                : -1.0 F (Oliver, present at event) / -0.645 R
%                           (FSG design doc). Three-way conflict on record.
%
% REQUIRES  hoosier_16x75_10_R20.m  (TTC-derived; camber derating lives
%           THERE and must not be applied again here)
% ------------------------------------------------------------------------

clear; clc; close all;

%% ---- parameters --------------------------------------------------------
P.m=306; P.L=1.530; P.tf=1.200; P.tr=1.190; P.g=9.80665;
P.ff=0.48;  P.ffBand=[0.43 0.53];
P.ms=P.m-2*8-2*13; P.hs=0.2993; P.muf=8; P.mur=13; P.hu=0.195;
P.RCf=0.0486; P.RCr=0.0533;
P.elasticFront=0.594;               % derived; reproduces 844/899 N.m/deg
P.rho=1.225;
P.aero=[3.07 0.490; 3.58 0.537];    % [SCz frontFrac], consistent pairs
P.aeroName={'handover','telemetry'};
P.IAf=1.0; P.IAr=0.645;             % deg magnitude
P.toe=[0 0.5 1.0];                  % deg PER SIDE
P.SC=0.585;                         % belt-to-track scaling, mid-range fit
P.V=11.5;                           % m/s, median steady-state speed

%% ---- 1. sub-limit balance ---------------------------------------------
ays=0.4:0.05:1.4;
fprintf('=== SUB-LIMIT BALANCE: alpha_f - alpha_r [deg], +ve understeer ===\n');
fprintf('%6s | %-22s | %-22s\n','','handover aero','telemetry aero');
fprintf('%6s %7s %7s %7s  %7s %7s %7s\n','a_y','toe0','toe.5','toe1','toe0','toe.5','toe1');
D=nan(numel(ays),3,2);
for i=1:numel(ays)
    for pA=1:2
        for j=1:3
            [af,ar]=axleAlphas(ays(i),P.V,P.ff,P.toe(j),pA,P);
            D(i,j,pA)=af-ar;
        end
    end
    if any(abs(ays(i)-[0.6 1.0 1.2])<1e-9)
        fprintf('%6.2f %+7.2f %+7.2f %+7.2f  %+7.2f %+7.2f %+7.2f\n', ...
            ays(i),D(i,1,1),D(i,2,1),D(i,3,1),D(i,1,2),D(i,2,2),D(i,3,2));
    end
end
fprintf(['-> weak OVERSTEER bias WITH THE AS-RUN TOE-OUT. At zero toe the same\n' ...
         '   model returns mild understeer, so the 1 deg front toe-out is the\n' ...
         '   MECHANISM, not a contributing factor.\n']);

Dlo=nan(size(ays)); Dhi=Dlo;
for i=1:numel(ays)
    v=[];
    for k=1:2, for pA=1:2
        [af,ar]=axleAlphas(ays(i),P.V,P.ffBand(k),0.5,pA,P); v(end+1)=af-ar; %#ok<AGROW>
    end, end
    Dlo(i)=min(v); Dhi(i)=max(v);
end

%% ---- 2. limiting axle across the uncertainty space --------------------
fprintf('\n=== LIMITING AXLE MAP (V=%.1f m/s, scaling %.3f, toe 0.5/side) ===\n',P.V,P.SC);
fprintf('%6s %11s %8s %11s %10s %8s %7s\n', ...
        'ff','aero','ay_max','front util','rear util','diff','first');
nRear=0; nTot=0; gaps=[];
for ff=[0.43 0.48 0.53]
    for pA=1:2
        [am,cf,cr,nf,nr]=limitState(P.V,ff,0.5,pA,P);
        uf=nf/cf; ur=nr/cr; nTot=nTot+1; nRear=nRear+(ur>=uf);
        gaps(end+1)=100*abs(uf-ur); %#ok<AGROW>
        fprintf('%6.2f %11s %8.3f %10.1f%% %9.1f%% %7.1f%% %7s\n', ...
            ff,P.aeroName{pA},am,100*uf,100*ur,100*abs(uf-ur), ...
            tern(ur>=uf,'REAR','FRONT'));
    end
end
fprintf(['-> rear-limited in %d of %d combinations; utilisation differs by\n' ...
         '   only %.1f-%.1f%%. The car is close to SIMULTANEOUSLY LIMITED and\n' ...
         '   the limiting axle is INDETERMINATE on this parameter set.\n'], ...
         nRear,nTot,min(gaps),max(gaps));

%% ---- 3. limit vs speed, both pairs, vs measured -----------------------
Vs=8:0.75:18;
ayA=arrayfun(@(v) limitState(v,P.ff,0.5,1,P), Vs);
ayB=arrayfun(@(v) limitState(v,P.ff,0.5,2,P), Vs);

mv=[]; mx=[];
CSV='IvanAxelEnduranceFSG25_2025Car_GenericTesting_a_3780.csv';
if isfile(CSV)
    T=readtable(CSV,'VariableNamingRule','preserve');
    Vv=T.("speed"); ay=T.("InlineAcc");
    ay=ay-mean(ay(Vv<0.3),'omitnan');
    ay=[nan;nan;ay(1:end-2)];                 % deskew +2 samples
    sus=movmin(abs(ay),10,'omitnan');         % sustained > 0.5 s
    for vlo=8:1.5:16.5
        k=Vv>=vlo & Vv<vlo+1.5;
        if nnz(k)<200, continue; end
        mv(end+1)=vlo+0.75; mx(end+1)=max(sus(k),[],'omitnan'); %#ok<AGROW>
    end
    % A measured sustained maximum is only a LIMIT if the car was actually
    % asked for one in that band. Bands where the demand never approached
    % capability (long straights, geometry-constrained hairpins) report a
    % low value that says nothing about the model. Flag them rather than
    % averaging them in.
    fprintf('\n=== PREDICTED vs MEASURED LIMIT ===\n');
    fprintf('  %8s %9s %11s %9s   %s\n','V m/s','model g','measured g','error','note');
    err=[]; atLimit=false(size(mv));
    for i=1:numel(mv)
        pr=limitState(mv(i),P.ff,0.5,2,P);
        e=100*(mx(i)/pr-1);
        atLimit(i)= e > -15;          % within 15% of predicted capability
        if atLimit(i), err(end+1)=e; note=''; %#ok<AGROW>
        else, note='limit not demanded in this band';
        end
        fprintf('  %8.1f %9.3f %11.3f %8.0f%%   %s\n',mv(i),pr,mx(i),e,note);
    end
    % ASYMMETRIC EVIDENCE: a sustained maximum is a LOWER BOUND on
    % capability. Bands where measurement falls below the model prove
    % nothing (the limit may simply not have been demanded). Only bands
    % where measurement EXCEEDS the model are binding - the car
    % demonstrably did it, so the model is too low by at least that much.
    binding = err(err>0);
    if isempty(binding)
        fprintf('  -> no band exceeds the model; measurement is a lower bound throughout.\n');
    else
        fprintf(['  -> BINDING RESULT: %d band(s) EXCEED the model, by up to %+.0f%%.\n' ...
                 '     The model under-predicts limit capability by at least that.\n' ...
                 '     Bands below the model are uninformative (limit not demanded)\n' ...
                 '     and are NOT evidence of agreement.\n'], numel(binding), max(binding));
    end
    fprintf(['     Where measurement EXCEEDS the model, the mid-range scaling\n' ...
             '     (0.585) understates limit capability - matching the measured\n' ...
             '     peak needs ~0.62-0.64. Candidate cause: combined-slip\n' ...
             '     contamination suppressing mid-range measured force. OPEN.\n']);
end

%% ---- figures -----------------------------------------------------------
GRN=[.26 .90 .58]; ORA=[1 .76 .20]; RED=[.94 .33 .31]; SFT=[.53 .60 .67];

figure('Color','k','Position',[60 60 1050 660]); hold on; grid on;
fill([ays fliplr(ays)],[Dlo fliplr(Dhi)],ORA,'FaceAlpha',0.12,'EdgeColor','none');
plot(ays,D(:,1,2),'-','Color',SFT,'LineWidth',1.4);
plot(ays,D(:,2,2),'-','Color',ORA,'LineWidth',2.6);
plot(ays,D(:,3,2),'-','Color',[1 .44 .26],'LineWidth',1.8);
plot(ays,D(:,2,1),'--','Color',ORA,'LineWidth',1.3);
yline(0,'Color',[.53 .53 .53]);
text(0.42,0.24,'UNDERSTEER (\alpha_f > \alpha_r)','Color',GRN,'FontWeight','bold');
text(0.42,-0.30,'OVERSTEER (\alpha_r > \alpha_f)','Color',RED,'FontWeight','bold');
xlabel('Lateral acceleration [g]'); ylabel('\alpha_f - \alpha_r  [deg]');
title({'TBRe25 sub-limit balance - tyre-model inversion, no steering channel', ...
 sprintf('V=%.1f m/s | camber -1.0F/-0.645R | scaling %.3f | elastic split 0.594F derived',P.V,P.SC)});
legend({'front fraction 43-53% \times both aero pairs','zero toe', ...
        'toe-out 0.5\circ/side (likely), telemetry aero','toe-out 1.0\circ/side', ...
        'toe-out 0.5\circ/side, handover aero'},'Location','southwest');
ylim([-0.8 0.45]); xlim([0.4 1.4]); dark(gca);

figure('Color','k','Position',[60 60 1050 660]); hold on; grid on;
plot(Vs,ayA,'--','Color',SFT,'LineWidth',1.8);
plot(Vs,ayB,'-','Color',GRN,'LineWidth',2.2);
if ~isempty(mv), plot(mv,mx,'o','Color',ORA,'MarkerFaceColor',ORA,'MarkerSize',8); end
xlabel('Speed [m/s]'); ylabel('Limit lateral acceleration [g]');
title({'TBRe25 predicted vs measured lateral limit', ...
       ['measured maxima are a LOWER BOUND on capability - only points ' ...
        'ABOVE the model are binding']});
legend({'model, handover aero','model, telemetry aero', ...
        'measured sustained |a_y| maxima (>0.5 s)'},'Location','northwest');
ylim([0.9 1.85]); dark(gca);
% NOTE: the lowest measured point (~0.95 g) is a band where the limit was
% never demanded. It is plotted, not hidden - dropping it would flatter the
% comparison.

%% ========================================================================
function [af,ar]=axleAlphas(ay,V,ff,toe,pA,P)
[Fzf0,Fzr0]=staticLoads(V,ff,pA,P);
[dFf,dFr]=loadTransfer(ay,ff,P);
af=invertAxle(P.m*ay*P.g*ff,     Fzf0,dFf,P.IAf,toe,P.SC);
ar=invertAxle(P.m*ay*P.g*(1-ff), Fzr0,dFr,P.IAr,0.0,P.SC);
end

function [Fzf0,Fzr0]=staticLoads(V,ff,pA,P)
% Per-tyre static + aero load. Aero from ONE consistent pair, computed per
% speed - downforce scales with v^2, so a fixed value is wrong away from
% the speed it was evaluated at.
SCz=P.aero(pA,1); aFF=P.aero(pA,2);
Dn=SCz*0.5*P.rho*V^2;
Fzf0=(P.m*P.g*ff+aFF*Dn)/2;
Fzr0=(P.m*P.g*(1-ff)+(1-aFF)*Dn)/2;
end

function [dFf,dFr]=loadTransfer(ay,ff,P)
% Elastic + geometric + unsprung, per tyre. Only the elastic part follows
% the roll-stiffness split; lumping all three into m*ay*h/t misallocates
% transfer between axles.
arm=P.hs-(P.RCf*(1-ff)+P.RCr*ff);
dFf=P.elasticFront    *P.ms*ay*P.g*arm/P.tf + P.ms*ff    *ay*P.g*P.RCf/P.tf + 2*P.muf*ay*P.g*P.hu/P.tf;
dFr=(1-P.elasticFront)*P.ms*ay*P.g*arm/P.tr + P.ms*(1-ff)*ay*P.g*P.RCr/P.tr + 2*P.mur*ay*P.g*P.hu/P.tr;
end

function alpha=invertAxle(Freq,Fz0,dF,IA,toe,SC)
if axleForce(14,Fz0,dF,IA,toe,SC)<Freq, alpha=NaN; return; end
lo=1e-5; hi=14;
for k=1:60
    mid=(lo+hi)/2;
    if axleForce(mid,Fz0,dF,IA,toe,SC)<Freq, lo=mid; else, hi=mid; end
end
alpha=(lo+hi)/2;
end

function F=axleForce(al,Fz0,dF,IA,toe,SC)
% Axle pair force. Toe splits the pair: loaded outer runs al+toe, inner
% al-toe. Magnitudes are summed because both tyres act in the same sense
% once the axle is generating net force; al-toe may go negative at small
% al, in which case that tyre opposes and its contribution is subtracted.
o=Fz0+dF; in=max(Fz0-dF,0);
F=abs(hoosier_16x75_10_R20(al+toe,o,SC,IA));
if in>1
    ai=al-toe;
    Fi=abs(hoosier_16x75_10_R20(ai,in,SC,IA));
    if ai>=0, F=F+Fi; else, F=F-Fi; end
end
end

function [aym,capF,capR,needF,needR]=limitState(V,ff,toe,pA,P)
lo=0.3; hi=3.5;
for k=1:48
    mid=(lo+hi)/2;
    [cf,cr,nf,nr]=state(mid,V,ff,toe,pA,P);
    if min(cf-nf,cr-nr)>0, lo=mid; else, hi=mid; end
end
aym=lo;
[capF,capR,needF,needR]=state(lo,V,ff,toe,pA,P);
end

function [capF,capR,needF,needR]=state(ay,V,ff,toe,pA,P)
[Fzf0,Fzr0]=staticLoads(V,ff,pA,P);
[dFf,dFr]=loadTransfer(ay,ff,P);
aa=linspace(0.2,14,240); capF=0; capR=0;
for al=aa
    capF=max(capF, axleForce(al,Fzf0,dFf,P.IAf,toe,P.SC));
    capR=max(capR, axleForce(al,Fzr0,dFr,P.IAr,0,  P.SC));
end
needF=P.m*ay*P.g*ff;
needR=P.m*ay*P.g*(1-ff);
end

function s=tern(c,a,b), if c, s=a; else, s=b; end, end
function dark(ax)
set(ax,'Color','k','XColor','w','YColor','w','GridColor',[.45 .45 .45],'GridAlpha',.4);
set(ax.Title,'Color','w');
lg=legend; if ~isempty(lg), set(lg,'TextColor','w','Color',[.12 .12 .12]); end
end
