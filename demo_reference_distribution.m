function out = demo_reference_distribution(showFigure)
%DEMO_REFERENCE_DISTRIBUTION Equal-amplitude theory curves in Fig. 7d,e.
% out = demo_reference_distribution(false) computes the results without plotting.
% G is the short-time mean change in Q = Rplus^2 + Rminus^2 due to noise,
% relative to noise-free evolution and divided by noise intensity and time.
if nargin < 1
    showFigure = true;
end
validateattributes(showFigure, {'logical'}, {'scalar'});
p = struct('J1',7,'J2',8,'K1',1.5,'K2',9,'Delta',1,'N',10000,'D',1);

% The 10,000-node Lorentzian quadrature uses each positive (v,w) and
% its three sign reflections. Reversing both frequencies leaves every
% required moment unchanged. Keep both signs of w to evaluate the two
% branch distributions separately at each prescribed amplitude pair.
v = p.Delta*tan(pi/2*radical_inverse(p.N/4,2));
w = p.Delta*tan(pi/2*radical_inverse(p.N/4,3));
omega = [v+w, v-w; v-w, v+w];
R = linspace(0.05,0.95,101)';
G = zeros(numel(R),2);
field = zeros(numel(R),2);
c2 = zeros(numel(R),2);
matchingError = zeros(size(R));
for k = 1:numel(R)
    target = R(k);
    % Initialize with the single-branch field. Both field strength and
    % second harmonic will be determined anew for the coupled branches.
    start = [log(4*target/(1-target^2))*[1;1]; target^2*[1;1]];
    [state,m,matchingError(k)] = match_amplitude(start,target,omega,p);
    field(k,:) = exp(state(1:2))';
    c2(k,:) = m.c2;
    for channel = 1:2
        G(k,channel) = order_response(m,channel);
    end
end

crossing = zeros(1,2);
for channel = 1:2
    j = find(G(1:end-1,channel)>=0 & G(2:end,channel)<0,1);
    crossing(channel) = interp1(G(j:j+1,channel),R(j:j+1),0);
end
out = struct('params',p,'R',R,'Gpair',G(:,1),'Ghigher',G(:,2), ...
    'crossing',crossing,'alpha',field,'c2',c2,'matchingError',matchingError);
if showFigure
    figure('Color','white','Name','Equal-amplitude reference distributions', ...
        'Position',[120 100 740 440]);
    ax = axes('FontName','Arial','FontSize',12,'TickDir','out','Box','off');
    hold(ax,'on');
    bounds = [min(G(:))*1.12, max(G(:))*1.35];
    patch(ax,crossing([1 2 2 1]),bounds([1 1 2 2]),[.92 .94 .91], ...
        'EdgeColor','none','HandleVisibility','off');
    text(ax,mean(crossing),bounds(2)-0.04*diff(bounds), ...
        'Opposite signs','HorizontalAlignment','center','Color',[.23 .32 .26], ...
        'FontName','Arial','FontSize',10);
    yline(ax,0,'-','Color',[.58 .59 .60],'HandleVisibility','off');
    colors = [33 102 172;196 91 40]/255;
    plot(ax,R,G(:,1),'-','Color',colors(1,:),'LineWidth',2.4, ...
        'DisplayName','Pairwise coupling noise');
    plot(ax,R,G(:,2),'--','Color',colors(2,:),'LineWidth',2.4, ...
        'DisplayName','Higher-order coupling noise');
    for channel = 1:2
        plot(ax,crossing(channel),0,'o','Color',colors(channel,:), ...
            'MarkerFaceColor','white','LineWidth',1.5,'HandleVisibility','off');
    end
    xlim(ax,[R(1) R(end)]); ylim(ax,bounds);
    xlabel(ax,'Initial amplitude   R_{+} = R_{-} = R');
    ylabel(ax,'Short-time mean response of Q');
    title(ax,'Noise response at equal branch amplitudes','FontWeight','normal');
    legend(ax,'Location','southwest','Box','off');
end
end

function g = order_response(m,channel)
% Eq. (15), evaluated before assuming exact equality of higher moments on
% the two branches. Eqs. (16)-(17) are its branch-symmetric forms.
opposite = [2 1];
sinSquared = (1-m.c2)/2;
if channel==1
    strength = m.R;
    dSum = strength.*sinSquared;
    dDifference = strength*m.sinCross;
else
    strength = m.R.*m.c2;
    % Derivatives of the three-body field along the sum and difference
    % of the two noise directions, including the changing second harmonic.
    dSum = 4*m.R.*strength.*m.cosSin2 + m.c2.*strength.*sinSquared;
    dDifference = 4*m.R(opposite).*strength.*m.cosSinCross(opposite) ...
        + m.c2(opposite).*strength*m.sinCross;
end
lambda = dSum.*sinSquared - 2*strength.^2.*m.cosSin2 ...
    + dDifference*m.sinCross - strength(opposite).^2.*m.cosOtherSin2 ...
    - strength(opposite).*strength.*m.cosOtherSin2(opposite);
bSum = strength.*sinSquared;
bDifference = strength(opposite)*m.sinCross;
g = sum(m.R.*lambda + bSum.^2 + bDifference.^2);
end

function q = radical_inverse(n,base)
% Radical-inverse quadrature in bases 2 and 3, starting with index 1.
index = (1:n)';
q = zeros(n,1);
weight = 1/base;
while any(index)
    q = q + weight*mod(index,base);
    index = floor(index/base);
    weight = weight/base;
end
end

function [state,m,err] = match_amplitude(start,target,omega,p)
% Adjust log(alpha) and c2 separately on each branch. Enforce
% <cos(phi)>=R and <cos(2*phi)>=c2 for the finite quadrature.
state = start;
[residual,m] = matching_residual(state,target,omega,p);
for k = 1:24
    if norm(residual,inf)<1e-9
        break
    end
    jac = zeros(4);
    for j = 1:4
        step = 2e-5;
        trial = state;
        trial(j) = trial(j)+step;
        jac(:,j) = (matching_residual(trial,target,omega,p)-residual)/step;
    end
    if rcond(jac)<1e-12
        break
    end
    correction = jac\residual;
    accepted = false;
    for j = 0:12
        trial = state - correction/(2^j);
        if any(trial(3:4)<-0.1 | trial(3:4)>1) || any(abs(trial(1:2))>12)
            continue
        end
        [r,mTrial] = matching_residual(trial,target,omega,p);
        if norm(r)<norm(residual)
            state = trial; residual = r; m = mTrial;
            accepted = true;
            break
        end
    end
    if ~accepted
        break
    end
end
err = norm(residual,inf);
if err>1e-7
    % Derivatives can be nonsmooth where an individual quadrature node locks.
    objective = @(x) norm(matching_residual(x,target,omega,p))^2;
    state = fminsearch(objective,state,optimset('Display','off', ...
        'TolX',1e-10,'TolFun',1e-16,'MaxIter',400));
    [residual,m] = matching_residual(state,target,omega,p);
    err = norm(residual,inf);
end
assert(err<2e-4,'Reference distribution could not be matched at R = %.4f.',target);
end

function [residual,m] = matching_residual(state,target,omega,p)
alpha = exp(state(1:2))';
c2 = state(3:4)';
% Ratio of cross-branch to same-branch fields in Eq. (S61).
beta = alpha.*((p.J1-p.K1)+(p.J2-p.K2)*c2) ./ ((p.J1+p.K1)+(p.J2+p.K2)*c2);
[s,locked] = conditional_sines(omega./alpha,beta([2 1])./alpha);
cosine = zeros(size(s));
cosine(locked) = sqrt(max(0,1-s(locked).^2));
harmonic = -s.^2;
harmonic(locked) = 1-2*s(locked).^2;
sinSquared = (1-harmonic)/2;
sinCross = prod(s,2);
m.c2 = mean(harmonic,1);
m.R = mean(cosine,1);
m.sinCross = mean(sinCross);
% These are conditional angular moments averaged over the SAME (v,w)
% nodes. A drifting angle gives zero to <cos(phi)*sin(phi)^2>.
m.cosSin2 = mean(cosine.*s.^2,1);
m.cosOtherSin2 = mean(cosine.*sinSquared(:,[2 1]),1);
m.cosSinCross = mean(cosine.*sinCross,1);
residual = [m.R-target, m.c2-c2]';
end

function [s,locked] = conditional_sines(f,ratio)
% f = Omega/alpha. The two columns are the + and - branches. Solve
% Eqs. (S62)-(S64) at fixed paired frequencies before taking averages.
bothValue = (f-ratio.*f(:,[2 1]))/(1-prod(ratio));
bothLocked = all(abs(bothValue)<1,2);
s = zeros(size(f));
s(bothLocked,:) = bothValue(bothLocked,:);
locked = [bothLocked,bothLocked];
assigned = bothLocked;
for j = 1:2
    % If branch j locks, eliminate its sine to solve the drifting
    % branch's quadratic. The smaller-magnitude root is physical.
    other = 3-j;
    v = f(:,other)-ratio(other)*f(:,j);
    discriminant = v.^2-1+2*prod(ratio);
    denom = v+sign(v).*sqrt(max(0,discriminant));
    driftSine = inf(size(v));
    regular = abs(denom)>1e-14;
    driftSine(regular) = 1./denom(regular);
    lockSine = f(:,j)-ratio(j)*driftSine;
    mixed = ~assigned & discriminant>0 & abs(driftSine)<1 ...
        & abs(lockSine)<1 & abs(f(:,other)-ratio(other)*lockSine)>1;
    s(mixed,j) = lockSine(mixed);
    s(mixed,other) = driftSine(mixed);
    locked(mixed,j) = true;
    assigned = assigned | mixed;
end

% Both branches drift for the remaining frequency pairs. Apply Eq. (S64)
% to the effective frequency on each branch until both means agree.
f = f(~assigned,:);
drift = zeros(size(f));
for k = 1:160
    u = f-ratio.*drift(:,[2 1]);
    next = u;
    free = abs(u)>1;
    next(free) = sign(u(free))./(abs(u(free))+sqrt(u(free).^2-1));
    change = max(abs(next-drift),[],'all');
    drift = (drift+next)/2;
    if isempty(change) || change<2e-11
        break
    end
end
s(~assigned,:) = drift;
end

