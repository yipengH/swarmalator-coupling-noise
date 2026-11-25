function wave = demo_phase_wave(showFigure)
%DEMO_PHASE_WAVE Local amplitude fluctuations of a continuum phase wave.
%   wave = demo_phase_wave() plots the response coefficients of Fig. 4d,e
%   and their prediction for the variance growth rate at the D of Fig. 4c.
%   wave = demo_phase_wave(false) returns the same values without plotting.
%
%   In Eq. (9), the angle x + theta is uniformly distributed and the
%   ordered branch has amplitude R = |mean(exp(i*(x - theta)))|.
%   The coefficients below follow from the Ott-Antonsen reduction.

if nargin == 0
    showFigure = true;
end
validateattributes(showFigure, {'logical'}, {'scalar'});

% Apply noise to one channel at a time with the same intensity D.
wave.D = 10;
wave.R = linspace(0,1,501)';
R = wave.R;
wave.a1 = 0.5 * R .* (1 - R.^2);
wave.a2 = 0.5 * R.^3 .* (1 - R.^2);

% Starting at R, the variance of a short amplitude increment divided
% by its duration approaches D*a^2 as the duration tends to zero.
wave.varianceRate1 = wave.D * wave.a1.^2;
wave.varianceRate2 = wave.D * wave.a2.^2;

% Check both coefficients from an angular distribution at R = 1/2.
% eta = x - theta is measured from the ordered branch's mean angle.
% Integrate the OA Poisson kernel on an equally spaced periodic grid.
check.R = 0.5;
check.numAngles = 128;
dEta = 2*pi / check.numAngles;
eta = (0:check.numAngles-1)' * dEta;
density = (1 - check.R^2) ./ (2*pi*(1 - 2*check.R*cos(eta) + check.R^2));
mass = density * dEta;
check.c1 = sum(mass .* cos(eta));
check.c2 = sum(mass .* cos(2*eta));
check.meanSinSquared = sum(mass .* sin(eta).^2);

% The pairwise and three-body fields are proportional to -sin(eta)
% with factors c1 and c1*c2. Projecting onto R = <cos(eta)> supplies
% a second factor of -sin(eta) before averaging over the angles.
check.a1 = check.c1 * check.meanSinSquared;
check.a2 = check.c1 * check.c2 * check.meanSinSquared;
atCheck = (R == check.R);
assert(abs(check.a1 - wave.a1(atCheck)) < 1e-12);
assert(abs(check.a2 - wave.a2(atCheck)) < 1e-12);
wave.check = check;

if showFigure
    pair = [33 102 172]/255;
    higher = [196 91 40]/255;
    f = figure('Name','Phase-wave fluctuations','Color','w','Position',[80 80 860 370]);
    panels = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
    ax = nexttile(panels);
    h1 = plot(ax,R,wave.a1,'-','Color',pair,'LineWidth',2.3);
    hold(ax,'on');
    h2 = plot(ax,R,wave.a2,'--','Color',higher,'LineWidth',2.3);
    % The maxima of the coefficients occur at different amplitudes.
    r1 = 1/sqrt(3); r2 = sqrt(3/5);
    plot(ax,r1,0.5*r1*(1-r1^2),'o','Color',pair,'MarkerFaceColor',pair);
    plot(ax,r2,0.5*r2^3*(1-r2^2),'o','Color',higher,'MarkerFaceColor',higher);
    xlim(ax,[0 1]); ylim(ax,[0 0.22]);
    xlabel(ax,'Phase-wave amplitude, R'); ylabel(ax,'Noise coefficient');
    title(ax,'Response of the ordered branch','FontWeight','normal');

    ax2 = nexttile(panels);
    plot(ax2,R,wave.varianceRate1,'-','Color',pair,'LineWidth',2.3);
    hold(ax2,'on');
    plot(ax2,R,wave.varianceRate2,'--','Color',higher,'LineWidth',2.3);
    xlim(ax2,[0 1]); ylim(ax2,[0 0.4]);
    xlabel(ax2,'Phase-wave amplitude, R'); ylabel(ax2,'Variance growth rate');
    title(ax2,'Short-time variance at D = 10','FontWeight','normal');
    set([ax ax2],'FontName','Arial','FontSize',11,'TickDir','out','Box','off');
    lg = legend(ax,[h1 h2],{'Pairwise coupling noise','Higher-order coupling noise'}, ...
        'Orientation','horizontal','Box','off');
    lg.Layout.Tile = 'south';
end
end
