function out = demo_first_passage(showFigure)
%DEMO_FIRST_PASSAGE Mean downward passage times in Fig. 6a,b (Eq. 11).
%   out = demo_first_passage(false) skips the figure. out.log10MFPT
%   contains the complete pairwise (column 1) and higher-order (column 2)
%   curves, including the times above the displayed observation window.
if nargin == 0
    showFigure = true;
end
validateattributes(showFigure, {'logical'}, {'scalar'});

p.J = [1, 5];
p.K = [6.5, 9];
p.Delta = 1;
p.epsilon = 0.01;                     % absorbing low-amplitude target
p.observationTime = 1000;
D = 10.^(-2:0.2:2)';
Jmean = (p.J+p.K)/2;
gamma = 2*p.Delta;                    % half-width of v - omega
q = sort(roots([-Jmean(2), Jmean(2)-Jmean(1), Jmean(1)-2*gamma]));
RH = sqrt(q(2));                      % starting amplitude

% The natural upper endpoint is 1. The last grid point is used only for
% quadrature. The code checks an upper bound on the mass omitted up to 1.
n = 12001;
low = exp(linspace(log(p.epsilon),log(RH),n))';
low([1 end]) = [p.epsilon; RH];
high = linspace(RH,1-1e-5,n)';
r = [low; high(2:end)];
h = diff(r);
weights = [h(1); h(1:n-2)+h(2:n-1); h(n-1)]/2; % outer trapezoid weights
log10MFPT = zeros(numel(D),2);
logTailFraction = zeros(numel(D),2);
for channel = 1:2
    a = 0.5*r.^(2*channel-1).*(1-r.^2);
    % V is an antiderivative of -F/a^2 in Eq. (10). Its constant
    % cancels from the inner integral through potential differences.
    v = drift_primitive(r,channel,Jmean,gamma);
    for j = 1:numel(D)
        U = log(a) + (2/D(j))*v;
        % At node y integrate exp[U(y)-U(z)] from z=y to the top.
        % Within each cell integrate the exponential of the linearly
        % interpolated potential exactly, including very steep cells.
        change = U(1:end-1)-U(2:end);
        magnitude = abs(change);
        logCell = log(h);
        k = magnitude>0;
        logCell(k) = logCell(k) + max(change(k),0) + ...
            log(-expm1(-magnitude(k))) - log(magnitude(k));
        logInner = -inf(size(r));
        for i = numel(h):-1:1
            rest = change(i) + logInner(i+1);
            logInner(i) = max(logCell(i),rest) + ...
                log1p(exp(-abs(logCell(i)-rest)));
        end
        outerWeights = log(weights) - 2*log(a(1:n));
        logMass = log_sum(outerWeights+logInner(1:n));
        log10MFPT(j,channel) = (log(2/D(j))+logMass)/log(10);

        % U increases between the last node and 1 for these parameters.
        % That interval contributes at most its width times the largest
        % remaining integrand. Compare this bound with the integral.
        logTailFraction(j,channel) = log(1-r(end)) - U(end) + ...
            log_sum(outerWeights+U(1:n)) - logMass;
    end
end
assert(all(logTailFraction(:) < log(1e-12)), 'Unresolved upper endpoint.');
out = struct('params',p,'D',D,'RH',RH,'log10MFPT',log10MFPT);
if ~showFigure
    return
end

figure('Color','w','Name','First passage from the phase wave', ...
    'Position',[100 100 740 460]);
ax = axes;
hold(ax,'on');
colors = [33 102 172; 196 91 40]/255;
names = {'Pairwise coupling','Higher-order coupling'};
styles = {'-o','-s'};
for channel = 1:2
    plot(ax,D,log10MFPT(:,channel),styles{channel}, ...
        'Color',colors(channel,:),'LineWidth',2.2,'MarkerSize',4, ...
        'MarkerFaceColor','w','DisplayName',names{channel});
end
yline(ax,log10(p.observationTime),'--','T_{obs} = 1000', ...
    'Color',[.28 .32 .35],'LineWidth',1.1, ...
    'LabelHorizontalAlignment','left','HandleVisibility','off');
set(ax,'XScale','log','FontName','Arial','FontSize',12,'Box','off', ...
    'TickDir','out','LineWidth',0.9,'YGrid','on','GridAlpha',.14);
xlim(ax,[D(1) D(end)]);
ylim(ax,[0 5]);                        % zoom on times around Tobs
xticks(ax,10.^(-2:2));
xlabel(ax,'Noise intensity  D');
ylabel(ax,'log_{10}(mean first-passage time)');
title(ax,'Time for a phase wave to reach low amplitude','FontSize',16);
subtitle(ax,sprintf('Starting amplitude %.3f,   target %.2f',RH,p.epsilon));
legend(ax,'Location','northeast','Box','off');
end

function v = drift_primitive(r,channel,J,gamma)
q = r.^2;
s = 1-q;
if channel == 1
    v = (4*gamma-2*J(1))*log(r) + (sum(J)-2*gamma)*log1p(-q) ...
        + 2*gamma./s;
else
    v = (12*gamma-2*sum(J))*log(r) + (sum(J)-6*gamma)*log1p(-q) ...
        + (sum(J)-4*gamma)./q + (J(1)-2*gamma)./(2*q.^2) + 2*gamma./s;
end
end

function value = log_sum(logValues)
shift = max(logValues);
value = shift + log(sum(exp(logValues-shift)));
end
