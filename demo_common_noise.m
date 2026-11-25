function out = demo_common_noise(makePlot)
%DEMO_COMMON_NOISE Recompute the single trajectories in Figs. S1d and S3d.
%   out = demo_common_noise() integrates and plots both trajectories.
%   out = demo_common_noise(false) returns their amplitudes without plotting.
%   One common noise drives all positions and a second drives all phases.
if nargin == 0
    makePlot = true;
end
validateattributes(makePlot, {'logical'}, {'scalar'});

p.N = 10000;
p.Delta = 1;                      % Lorentzian half-width for v and omega
p.J1 = 2; p.J2 = 5;               % Position coupling strengths
p.K1 = 7; p.K2 = 9;               % Phase coupling strengths
p.D = 10;                         % Intensity in the fluctuating channel
p.dt = 0.01;
p.T = 120;
p.recordDt = 0.1;
p.seed = 2027102610;              % Shared by the two published trajectories

folder = fileparts(mfilename('fullpath'));
initial = load(fullfile(folder,'inputs','asynchronous_particles.mat'));
assert(all([numel(initial.x),numel(initial.theta), ...
    numel(initial.v),numel(initial.omega)] == p.N));

% Resetting the seed for each run pairs the realizations across channels.
originalRng = rng;
restoreRng = onCleanup(@() rng(originalRng));
out.params = p;
out.t = (0:round(p.T/p.recordDt))' * round(p.recordDt/p.dt) * p.dt;
out.pairwise = simulate(initial,p,1);
out.higherOrder = simulate(initial,p,2);
if makePlot
    figure('Name','Single particle-system realizations','Color','w', ...
        'Units','pixels','Position',[80 80 850 560]);
    layout = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
    ax = nexttile(layout);
    plot_branches(ax,out.t,out.pairwise,'Pairwise coupling noise');
    labels = legend(ax,{'R_{+}','R_{-}'},'Orientation','horizontal','Box','off');
    labels.Layout.Tile = 'south';
    ax = nexttile(layout);
    plot_branches(ax,out.t,out.higherOrder,'Higher-order coupling noise');
    xlabel(ax,'Time, t');
end
end

function series = simulate(initial,p,channel)
rng(p.seed,'twister');
x = initial.x;
theta = initial.theta;
recordEvery = round(p.recordDt/p.dt);
nSteps = round(p.T/p.dt);
series.Rplus = zeros(nSteps/recordEvery+1,1);
series.Rminus = zeros(size(series.Rplus));
[series.Rplus(1),series.Rminus(1)] = amplitudes(x,theta);
sigma = sqrt(2*p.dt*p.D);
for step = 1:nSteps
    [dx0,dt0,fx0,ft0] = coefficients(x,theta,initial,p,channel);
    noiseX = sigma*randn;
    noiseTheta = sigma*randn;
    xPredict = x + p.dt*dx0 + noiseX*fx0;
    thetaPredict = theta + p.dt*dt0 + noiseTheta*ft0;
    [dx1,dt1,fx1,ft1] = coefficients(xPredict,thetaPredict,initial,p,channel);

    % Stratonovich Heun corrector, using the same two common increments.
    x = x + 0.5*p.dt*(dx0+dx1) + 0.5*noiseX*(fx0+fx1);
    theta = theta + 0.5*p.dt*(dt0+dt1) + 0.5*noiseTheta*(ft0+ft1);
    if mod(step,100) == 0
        x = mod(x+pi,2*pi)-pi;
        theta = mod(theta+pi,2*pi)-pi;
    end
    if mod(step,recordEvery) == 0
        k = 1+step/recordEvery;
        [series.Rplus(k),series.Rminus(k)] = amplitudes(x,theta);
    end
end
end

function [driftX,driftTheta,noiseX,noiseTheta] = coefficients(x,theta,initial,p,channel)
[pairPlus,threePlus] = branch_fields(x+theta);
[pairMinus,threeMinus] = branch_fields(x-theta);
pairX = 0.5*(pairPlus+pairMinus);
pairTheta = 0.5*(pairPlus-pairMinus);
threeX = 0.5*(threePlus+threeMinus);
threeTheta = 0.5*(threePlus-threeMinus);
driftX = initial.v + p.J1*pairX + p.J2*threeX;
driftTheta = initial.omega + p.K1*pairTheta + p.K2*threeTheta;
if channel == 1
    noiseX = pairX;
    noiseTheta = pairTheta;
else
    noiseX = threeX;
    noiseTheta = threeTheta;
end
end

function [pair,three] = branch_fields(angle)
% Exact all-to-all pair and three-body sums using two complex moments.
z = exp(1i*angle);
Z1 = mean(z);
Z2 = mean(z.^2);
pair = imag(Z1.*conj(z));
three = imag(Z2.*conj(Z1).*conj(z));
end

function [rp,rm] = amplitudes(x,theta)
rp = abs(mean(exp(1i*(x+theta))));
rm = abs(mean(exp(1i*(x-theta))));
end

function plot_branches(ax,t,series,label)
% As in Figs. S1 and S3, colour and line style distinguish branches.
plot(ax,t,series.Rplus,'-','Color',[0.77 0.29 0.19],'LineWidth',1.25);
hold(ax,'on');
plot(ax,t,series.Rminus,'--','Color',[0.13 0.40 0.67],'LineWidth',1.4);
xlim(ax,[0 t(end)]); ylim(ax,[0 1]);
xticks(ax,0:20:120); yticks(ax,0:0.25:1);
ylabel(ax,'Branch amplitude');
title(ax,label,'FontWeight','normal');
set(ax,'FontName','Arial','FontSize',12,'TickDir','out','Box','off', ...
    'YGrid','on','GridAlpha',0.12,'Layer','top');
end
