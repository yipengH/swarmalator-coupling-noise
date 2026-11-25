function out = demo_probability_potential(showFigure)
%DEMO_PROBABILITY_POTENTIAL Local phase-wave potentials in Fig. 5b,e.
%   out = demo_probability_potential(false) skips the figure.
%   U(R) = -log[p(R)/p(Ru)] for the interior zero-current solution.
%   Columns of out.U are the pairwise and higher-order channels.
if nargin == 0
    showFigure = true;
end
validateattributes(showFigure, {'logical'}, {'scalar'});

p.J = [1, 5];                         % position couplings J1, J2
p.K = [6.5, 9];                       % phase couplings K1, K2
p.Delta = 1;
p.D = 1;                              % active channel; the other has D = 0
Jmean = (p.J + p.K)/2;
drift = @(r) -2*p.Delta*r + ...
    0.5*(Jmean(1) + Jmean(2)*r.^2).*r.*(1-r.^2);
q = sort(roots([-Jmean(2), Jmean(2)-Jmean(1), Jmean(1)-4*p.Delta]));
Ru = sqrt(q(1));                      % noise-free unstable amplitude
RH = sqrt(q(2));                      % noise-free stable high amplitude
R = unique([linspace(0.01,0.999,20001)'; Ru; RH]);
atRu = find(R == Ru);
U = zeros(numel(R),2);
peak = zeros(1,2);
valley = zeros(1,2);
for channel = 1:2
    a = 0.5*R.^(2*channel-1).*(1-R.^2);
    % Integrate away from Ru in each direction to keep the zero exact.
    integrand = drift(R)./a.^2;
    left = atRu:-1:1;
    right = atRu:numel(R);
    primitive = zeros(size(R));
    primitive(left) = cumtrapz(R(left), integrand(left));
    primitive(right) = cumtrapz(R(right), integrand(right));
    U(:,channel) = log(a/a(atRu)) - (2/p.D)*primitive;

    % At D = 1 each channel has one peak followed by a high-R minimum.
    slope = diff(U(:,channel));
    peak(channel) = find(slope(1:end-1)>0 & slope(2:end)<0,1) + 1;
    valley(channel) = find(slope(1:end-1)<0 & slope(2:end)>0,1) + 1;
end
out = struct('params',p,'R',R,'Ru',Ru,'RH',RH,'U',U, ...
    'peakR',R(peak),'valleyR',R(valley), ...
    'deltaU',[U(peak(1),1)-U(valley(1),1), U(peak(2),2)-U(valley(2),2)]);
if ~showFigure
    return
end

figure('Color','w','Name','Phase-wave probability potential', ...
    'Position',[100 100 900 445]);
layout = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
title(layout,'Probability potential of the phase wave', ...
    'FontName','Arial','FontSize',17);
names = {'Pairwise coupling','Higher-order coupling'};
colors = [33 102 172; 196 91 40]/255;
limits = [-0.25 0.45; -8.5 1.5];    % show the two extrema on each curve
for channel = 1:2
    ax = nexttile(layout);
    hold(ax,'on');
    yline(ax,0,':','Color',[.62 .65 .67]);
    plot(ax,R,U(:,channel),'Color',colors(channel,:),'LineWidth',2.4);
    for index = [peak(channel), valley(channel)]
        plot(ax,R(index),U(index,channel),'o','MarkerSize',7, ...
            'MarkerFaceColor','w','MarkerEdgeColor',colors(channel,:), ...
            'LineWidth',1.7);
    end
    xlim(ax,[.25 .8]);
    ylim(ax,limits(channel,:));
    xlabel(ax,'Phase-wave amplitude  R');
    ylabel(ax,'Relative potential  U(R)');
    title(ax,names{channel},'FontSize',13);
    subtitle(ax,sprintf('D = 1,   \\DeltaU = %.3g',out.deltaU(channel)));
    set(ax,'FontName','Arial','FontSize',11,'Box','off','TickDir','out', ...
        'LineWidth',0.9,'YGrid','on','GridAlpha',.12);
end
end
