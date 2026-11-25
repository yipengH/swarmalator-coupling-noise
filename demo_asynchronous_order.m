function out = demo_asynchronous_order(showFigure)
%DEMO_ASYNCHRONOUS_ORDER Stationary small-amplitude theory in Fig. 3c,d.
%   out = demo_asynchronous_order(false) returns the curves without plotting.
if nargin == 0
    showFigure = true;
end
validateattributes(showFigure, {'logical'}, {'scalar'});

% Asynchronous row in Table 1. At leading order only J1 + K1 enters.
p.N = 10000;
p.Delta = 1;
p.J1 = 7; p.J2 = 8;
p.K1 = -1; p.K2 = 9;
p.densityD = 1;                     % Bottom row of Fig. 3d
gap = 8*p.Delta - p.J1 - p.K1;     % Positive on the stable side
out.params = p;
out.D = unique([logspace(-2,0,401), linspace(0.01,0.02,201)])';

% Eq. (8). Evaluate the gamma ratio in log space at weak noise.
shape = gap ./ out.D;
out.meanPairwise = exp(0.5*log(8*p.Delta ./ (p.N*out.D)) ...
    + gammaln(1.5) + gammaln(shape - 0.5) - gammaln(shape));
out.meanHigherOrder = sqrt(2*pi*p.Delta/(p.N*gap))*ones(size(out.D));

% Eq. (7). The two branches have the same marginal density at this order.
out.R = linspace(0,0.14,1601)';
r = out.R;
prefactor = (p.N*gap/(4*p.Delta))*r;
out.pdfPairwise = prefactor .* exp(-(1 + gap/p.densityD) ...
    *log1p(p.N*p.densityD*r.^2/(8*p.Delta)));
out.pdfHigherOrder = prefactor .* exp(-p.N*gap*r.^2/(8*p.Delta));

if ~showFigure
    return
end
pair = [33 102 172]/255;
higher = [196 91 40]/255;
f = figure('Name','Asynchronous order','Color','w','Position',[80 80 860 370]);
panels = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(panels);
h1 = semilogx(ax,out.D,out.meanPairwise,'-','Color',pair,'LineWidth',2.3);
hold(ax,'on');
h2 = semilogx(ax,out.D,out.meanHigherOrder,'--','Color',higher,'LineWidth',2.3);
xlim(ax,[0.01 1]); ylim(ax,[0.017 0.023]);
xticks(ax,[0.01 0.1 1]); xticklabels(ax,{'0.01','0.1','1'});
ytickformat(ax,'%.3f');
xlabel(ax,'Noise intensity, D'); ylabel(ax,'Mean amplitude');
title(ax,'Mean order near the asynchronous state','FontWeight','normal');

ax2 = nexttile(panels);
% R=0 has exactly zero density and lies outside a logarithmic y axis.
positive = r > 0;
semilogy(ax2,r(positive),out.pdfPairwise(positive),'-','Color',pair,'LineWidth',2.3);
hold(ax2,'on');
semilogy(ax2,r(positive),out.pdfHigherOrder(positive),'--','Color',higher,'LineWidth',2.3);
xlim(ax2,[0 0.12]); ylim(ax2,[0.03 90]); xticks(ax2,0:0.03:0.12);
xlabel(ax2,'Single-branch amplitude, R'); ylabel(ax2,'Probability density');
title(ax2,'Amplitude distribution at D = 1','FontWeight','normal');
set([ax ax2],'FontName','Arial','FontSize',11,'TickDir','out','Box','off');
lg = legend(ax,[h1 h2],{'Pairwise coupling noise','Higher-order coupling noise'}, ...
    'Orientation','horizontal','Box','off');
lg.Layout.Tile = 'south';
end
