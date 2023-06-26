%% Oscillatory Bioreactor Dynamics
% Alexis Dubs adubs@mit.edu
% 6/3/2023


% timespan for simulations
tspan = [1, 400];

% parameters
p.mumax = 0.3;
p.ks = 1.75;
p.vmin = 0.1;
p.a = 0.03;
p.rho = 0.75;

% function handles for yield and growth functions
v = @yield;
mu = @growth;

D = 0.25; % dilution rate
sf = 200; % feed substrate concentration
x0 = [sf,1]; %initial substrate and biomass concentrations

% equilibrium concentrations
s_eq = D * p.ks / (p.mumax - D); % equilibrium value of s
x_eq = (p.vmin + p.a*s_eq)^p.rho * (sf-s_eq); %eq value of x
s_wo = sf; %washout value of s
x_wo = 0; %washout value of x

%% Model OG model
% run model
[t, x] = ode15s(@(t, x) model(t, x, D, sf, mu, v, p), tspan, x0);

% calculate things throughout simulation time
phis = zeros(length(t),1);
mus = zeros(length(t),1);
vs = zeros(length(t),1);
phiprimes = zeros(length(t),1);
for i = 1:length(t)
    mus(i) = mu(x(i,1), p);
    vs(i) = v(x(i,1), p);
    phis(i) = mu(x(i,1), p)/v(x(i,1), p);
    phiprimes(i) = dphids(x(i,1), p);
end

% plot time series of concentration
figure
plot(t, x, 'LineWidth', 1)
hold on
% plot equibliriums
plot(t, s_eq*ones(length(t),1), 'LineWidth', 1, 'LineStyle','--', 'Color','#0072BD')
plot(t, x_eq*ones(length(t),1), 'LineWidth', 1, 'LineStyle','--','Color', "#D95319")
legend({'s', 'x', 's_{eq}', 'x_{eq}'}, 'Location', 'northwest')
xlabel('Time (hr)')
ylabel('Concentration (g/L)')
exportgraphics(gcf, 'ogmodel.png', 'Resolution', 300)

% plot time series of phi, mu, v, and dphi/ds
figure
plot(t, mus, 'LineWidth', 1)
hold on
plot(t, vs, 'LineWidth', 1)
plot(t, phiprimes, 'LineWidth', 1)
plot(t, phis, 'LineWidth', 1)
xlabel('Time (hr)')
legend({'\mu', 'v', 'd\phi/ds', '\phi'}, 'NumColumns', 2)
exportgraphics(gcf, 'ogmodel_phi.png', 'Resolution', 300)

% plot fancy phase portrait
xs = linspace(1, 350, 15);
ss = linspace(1, 200, 15);
[S, X] = meshgrid(ss, xs);
dsdt = zeros(size(S));
dxdt = zeros(size(S));
t0 = 0;
% calculate phase field
for i = 1:numel(S)
    yprime = model(t0,[S(i), X(i)], D, sf, mu, v, p);
    dsdt(i) = yprime(1);
    dxdt(i) = yprime(2);
end
figure
% plot phase field
q = quiver(S,X,dsdt, dxdt);
q.AutoScaleFactor = 2;
q.LineWidth = 1;
xlim([0,200]);
ylim([0,350]);
hold on
% plot system trajectory
plot(x(:,1), x(:,2), 'LineWidth', 1);
% plot equilibrium points
plot([s_eq, s_wo], [x_eq, x_wo], 'LineWidth', 1, 'Marker','*', 'Color', "#7E2F8E", 'MarkerSize',10, 'LineStyle', 'none');
xlabel('Substrate Concentration (g/L)')
ylabel('Biomass Concentration (g/L)')
legend({'Phase Plane', 'System Trajectory', 'Equilibrium Points'})
exportgraphics(gcf, 'phase_portrait.png', 'Resolution', 300)

%% Plot s vs phi to show bend
ss = linspace(0,50,500);
phis = zeros(length(ss),1);
for i = 1:length(ss)
    phis(i) = mu(ss(i), p)/v(ss(i), p);
end
figure
plot(ss, phis, 'LineWidth', 1)
xlabel('s')
ylabel('\phi(s)')
exportgraphics(gcf, 'phi.png', 'Resolution', 300)

%% Model OG model with control
Kc = 10; % controller gain
[t, x] = ode45(@(t, x) model_control(t, x, D, sf, mu, v, s_eq, Kc, p), tspan, x0);

% calculate sf during simulation
sfs = zeros(length(t),1);
for i = 1:length(t)
    sfs(i) = max(0,sf + Kc*(s_eq-x(i,1)));
end

% plot concentration time series
figure
plot(t, x, 'LineWidth', 1)
hold on
plot(t,sfs, 'LineWidth', 1)
plot(t, s_ss*ones(length(t),1), 'LineWidth', 1, 'LineStyle','--', 'Color','#0072BD')
plot(t, sf*ones(length(t),1), 'LineWidth', 1, 'LineStyle','--','Color','#EDB120')
xlabel('Time (hr)')
ylabel('Concentration (g/L)')
ylim([0,sf + 1])
legend({'s', 'x', 's_f', 's_{eq}', 's_{nom}'}, 'Location', 'northeast')
exportgraphics(gcf, 'ogmodel_control.png', 'Resolution', 300)


%% Change s_f
Kc = 10;
% run simulation
[t, x] = ode15s(@(t, x) model_control(t, x, D, sf, mu, v, s_eq, Kc, p), tspan, x0);
% calculate sf during simulation
sfs = zeros(length(t),1);
for i = 1:(length(t))
    sfs(i) = max(0,sf + Kc*(s_ss-x(i,1)));
end

% change sf/snom and run simulation
sf2 = 300;
[t2, x2] = ode15s(@(t, x) model_control(t, x, D, sf2, mu, v, s_eq, Kc, p), tspan, x(end,:));
% calculate sf during simulation
sfs2 = zeros(length(t2),1);
for i = 1:(length(t2))
    sfs2(i) = max(0,sf2 + 10*(s_ss-x2(i,1)));
end

% append things from both simulations
snoms = [sf*ones(length(t),1);sf2*ones(length(t2),1)];
t = [t;t(end)+t2];
x = [x;x2];
sfs = [sfs; sfs2];

% plot everything
figure
plot(t, x, 'LineWidth', 1)
hold on
plot(t,sfs, 'LineWidth', 1)
plot(t, s_ss*ones(length(t),1), 'LineWidth', 1, 'LineStyle','--', 'Color','#0072BD')
plot(t, snoms, 'LineWidth', 1, 'LineStyle','--','Color','#EDB120')
xlabel('Time (hr)')
ylabel('Concentration (g/L)')
legend({'s', 'x', 's_f', 's_{eq}', 's_{nom}'}, 'Location', 'northwest')
exportgraphics(gcf, 'change_sf.png', 'Resolution', 300)

%% Perfusion with different alphas
alphas = [0.5, 0.8, 1];
xs = cell(length(alphas),1);
ts = cell(length(alphas),1);

% run simulation and plot phase diagram for each alpha
figure
for i = 1: length(alphas)
    [ts{i}, xs{i}] = ode45(@(t, x) model_perf(t, x, D, sf, mu, v, alphas(i), p), tspan, x0);
    plot3(alphas(i)*ones(length(ts{i}),1),xs{i}(:,1),xs{i}(:,2))
    hold on
end
xlabel('\alpha')
ylabel('s (g/L)')
zlabel('x (g/L)')
exportgraphics(gcf, 'perf_phase.png', 'Resolution', 300)

% plot time series of substrate for each alpha
figure
for i = 1: length(alphas)
    plot(ts{i},xs{i}(:,1), 'LineWidth', 1, 'DisplayName', strcat('\alpha= ', num2str(alphas(i))))
    hold on
end
xlabel('Time (hr)')
ylabel('Substrate Concentration (g/L)')
legend
exportgraphics(gcf, 'perf_alphas_s.png', 'Resolution', 300)

% plot time series of biomass for each alpha
figure
for i = 1: length(alphas)
    plot(ts{i},xs{i}(:,2), 'LineWidth', 1, 'DisplayName', strcat('\alpha= ', num2str(alphas(i))))
    hold on
end
xlabel('Time (hr)')
ylabel('Biomass Concentration (g/L)')
legend
exportgraphics(gcf, 'perf_alphas_x.png', 'Resolution', 300)

%% Plot stability criteria
num = 1000;
ds = linspace(0, 0.48, num);
sfs = [1, 200];
num2 = length(sfs);
alpha0 = 0.5;
alphas = zeros(num,num2);

% calculate alpha value for stability boundary
for i = 1:num
    for j = 1:num2
        opts = optimoptions('fsolve', 'FunctionTolerance', eps, 'Display', 'off');
        result = fsolve(@(alpha) criteria(ds(i), alpha, sfs(j), p), alpha0, opts);
        if ~isreal(result) || result > 1
            alphas(i,j) = NaN;
        else
            alphas(i,j) = result;
        end
    end
end

% plot
figure
for i = 1:num2
    plot(ds,alphas(:,i), 'LineWidth', 1, 'DisplayName', strcat('s_f= ', num2str(sfs(i)), ' g/L'))
    hold on
end
xlabel('D (1/hr)')
ylabel('\alpha')
legend
exportgraphics(gcf, 'stab_crit.png', 'Resolution', 300)

%% Functions
function out = dphids(s, p)
% derivative of mu/v wrt s
    mumax = p.mumax;
    ks = p.ks;
    vmin = p.vmin;
    a = p.a;
    rho = p.rho;

    out = mumax*(vmin+a*s)^(-rho-1)*(ks*(vmin-a*(rho-1)*s)-a*rho*s^2)/(ks+s)^2;
end

function output = criteria(D, alpha, sf, p)
% value of alpha on boundary of stable and unstable system
    mumax = p.mumax;
    ks = p.ks;
    vmin = p.vmin;
    a = p.a;
    rho = p.rho;

    seq = alpha*D*ks/(mumax - alpha*D);
    xeq = (vmin + a*seq)^rho/alpha*(sf-seq);
    dphids = mumax*(vmin+a*seq)^(-rho-1)*(ks*(vmin-a*(rho-1)*seq)-a*rho*seq^2)/(ks+seq)^2;

    output = D + dphids*xeq;
end

function dxdt = model_perf(t, x, D, sf, mu, v, alpha, p)
    dxdt = zeros(2,1);
    s = x(1);
    X = x(2);
    dxdt(1) = D * (sf-s) - mu(s, p)/v(s, p) * X;
    dxdt(2) = (mu(s, p) - alpha*D) * X;
end

function dxdt = model(t, x, D, sf, mu, v, p)
    dxdt = zeros(2,1);
    s = x(1);
    X = x(2);
    dxdt(1) = D * (sf-s) - mu(s, p)/v(s, p) * X;
    dxdt(2) = (mu(s, p) - D) * X;
end

function dxdt = model_control(t, x, D, sf, mu, v, s_ss, Kc, p)
    dxdt = zeros(2,1);
    s = x(1);
    X = x(2);
    sf = max(0, sf + Kc*(s_ss-s));
    dxdt(1) = D * (sf-s) - mu(s, p)/v(s, p) * X;
    dxdt(2) = (mu(s, p) - D) * X;
end

function output = yield(s, p)
    output = (p.vmin + p.a * s)^p.rho;
end

function output = growth(s, p)
    output = p.mumax*s/(p.ks + s);
end