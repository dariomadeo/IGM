clear all
close all

%% Download COVID19 data from italian Protezione Civile github account

t1 = datetime(2020,03,24,8,0,0); % First day 
t2 = datetime(2021,9,14,8,0,0); % Last day

tInterval = t1:t2;
NInterval = numel(tInterval); % Number of days

covidTable = cell(NInterval, 1);

for i=1:NInterval

    date_str = datestr(tInterval(i), 'yyyymmdd');
    fname = downloadItalianCovidData(date_str, './data');
    covidTable{i} = readtable(fname, 'ReadVariableNames',true); 

end

% Extract the relevant data
hospitalized    = zeros(NInterval, 1);
quarantined     = zeros(NInterval, 1);
recovered       = zeros(NInterval, 1);
deceased        = zeros(NInterval, 1);

for i=1:NInterval

    ct = covidTable{i};

    for row=1:21 % Italian regions
        hospitalized(i) = hospitalized(i)   + ct.totale_ospedalizzati(row);
        quarantined(i)  = quarantined(i)    + ct.isolamento_domiciliare(row);
        recovered(i)    = recovered(i)      + ct.dimessi_guariti(row);
        deceased(i)     = deceased(i)       + ct.deceduti(row);

    end
end

Idata = hospitalized + quarantined;
Rdata = recovered;
Ddata = deceased;

Tw = 30; % Length of a time window
W = floor(numel(Idata)/Tw); % Number of time windows

%% General parameters
M = 60244639; % Italian population size

%% Parameters setup (IGM)
par_obj_IGM = setParameter(); % Create an empty parameter object

par_obj_IGM = setParameter(par_obj_IGM, 'S0', 'tinvar', 59*1e6, 0.98*M, M);
par_obj_IGM = setParameter(par_obj_IGM, 'x0', 'tinvar', 0.41, 0.1, 0.9);

par_obj_IGM = setParameter(par_obj_IGM, 'beta0', 'tinvar', 0.5, 0.1, 0.7);
par_obj_IGM = setParameter(par_obj_IGM, 'gamma', 'const', 1/14);

par_obj_IGM = setParameter(par_obj_IGM, 'sigmaC0', 'tinvar', 4.22, 1e-2, 10);
par_obj_IGM = setParameter(par_obj_IGM, 'sigmaN0', 'tinvar', -4.2, -10, -1e-2);

par_obj_IGM = setParameter(par_obj_IGM, 'alpha', 'tvar', 0.03, 1e-4, 0.2833);
par_obj_IGM = setParameter(par_obj_IGM, 'tau', 'const', 0.166);
par_obj_IGM = setParameter(par_obj_IGM, 'omega', 'tvar', 0.019, 1e-4, 1/7);
par_obj_IGM = setParameter(par_obj_IGM, 'zeta', 'tvar', 0.01, 1e-4, 0.3);

par_obj_IGM = setParameter(par_obj_IGM, 'e', 'tvar', 0.75, 0.05, 0.95);
par_obj_IGM = setParameter(par_obj_IGM, 'a', 'tvar', 87370, 1, 1e6);

[theta0_IGM, lb_IGM, ub_IGM] = getParameter(par_obj_IGM, W); % % Build up ic, lb and ub vectors for the optimization problem

%% Parameters setup (IM)
par_obj_IM = setParameter(); % Create an empty parameter object

par_obj_IM = setParameter(par_obj_IM, 'S0', 'tinvar', 59*1e6, 0.98*M, M);
par_obj_IM = setParameter(par_obj_IM, 'x0', 'const', 0); % n.a.

par_obj_IM = setParameter(par_obj_IM, 'beta0', 'tvar', 0.5, 0.1, 0.7); % n.a.
par_obj_IM = setParameter(par_obj_IM, 'gamma', 'const', 1/14); 

par_obj_IM = setParameter(par_obj_IM, 'sigmaC0', 'const', 0); % n.a.
par_obj_IM = setParameter(par_obj_IM, 'sigmaN0', 'const', 0); % n.a.

par_obj_IM = setParameter(par_obj_IM, 'alpha', 'tvar', 0.03, 1e-4, 0.2833);
par_obj_IM = setParameter(par_obj_IM, 'tau', 'const', 0.166);
par_obj_IM = setParameter(par_obj_IM, 'omega', 'tvar', 0.019, 1e-4, 1/7);
par_obj_IM = setParameter(par_obj_IM, 'zeta', 'tvar', 0.01, 1e-4, 0.3);

par_obj_IM = setParameter(par_obj_IM, 'e', 'const', 0); % n.a.
par_obj_IM = setParameter(par_obj_IM, 'a', 'const', 0); % n.a.

[theta0_IM, lb_IM, ub_IM] = getParameter(par_obj_IM, W); % Build up ic, lb and ub vectors for the optimization problem

%% Estimate
thetaOptim_IGM = estimate(W, Tw, Idata, Rdata, Ddata, lb_IGM, ub_IGM, theta0_IGM, par_obj_IGM, M);

thetaOptim_IM = estimate(W, Tw, Idata, Rdata, Ddata, lb_IM, ub_IM, theta0_IM, par_obj_IM, M);

%% Simulate
tBreaks = zeros(W, 1);
tBreaks(1) = Tw-1;
for i=2:W
    tBreaks(i) = tBreaks(i-1) + Tw;
end

[S_IGM, U_IGM, I_IGM, R_IGM, D_IGM, x_IGM] = simode(thetaOptim_IGM, ...
                                                    Idata, ...
                                                    Rdata, ...
                                                    Ddata, ...
                                                    M, ...
                                                    tBreaks, ...
                                                    par_obj_IGM);

[S_IM, U_IM, I_IM, R_IM, D_IM, x_IM] = simode(thetaOptim_IM, ...
                                              Idata, ...
                                              Rdata, ...
                                              Ddata, ...
                                              M, ...
                                              tBreaks, ...
                                              par_obj_IM);


%% Plot

colorlist = lines(7);

figure
subplot(3,1,1)
hold on
plot(0:tBreaks(end), Idata(1:tBreaks(end)+1), 'o', 'LineWidth', 2, 'Color', colorlist(3,:), 'MarkerSize', 5);
plot(0:tBreaks(end), I_IM(1:tBreaks(end)+1), 'LineWidth', 2, 'Color', colorlist(2,:));
plot(0:tBreaks(end), I_IGM(1:tBreaks(end)+1), 'LineWidth', 2, 'Color', colorlist(1,:));

xlim([0 539]) 
legend('Real data', 'Sim. IM', 'Sim. IGM', 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'northwest')
ylabel('$I(t)$', 'Interpreter', 'latex', 'FontSize', 12)
set(gca, 'TickLabelInterpreter', 'latex', 'fontsize', 12)
xlabel('$t$ (days)' , 'Interpreter', 'latex', 'FontSize', 12)

subplot(3,1,2)
hold on
plot(0:tBreaks(end), Rdata(1:tBreaks(end)+1), 'o', 'LineWidth', 2, 'Color', colorlist(3,:), 'MarkerSize', 5);
plot(0:tBreaks(end), R_IM(1:tBreaks(end)+1), 'LineWidth', 2, 'Color', colorlist(2,:));
plot(0:tBreaks(end), R_IGM(1:tBreaks(end)+1), 'LineWidth', 2, 'Color', colorlist(1,:));

xlim([0 539]) 
legend('Real data', 'Sim. IM', 'Sim. IGM', 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'northwest')
ylabel('$R(t)$', 'Interpreter', 'latex', 'FontSize', 12)
set(gca, 'TickLabelInterpreter', 'latex', 'fontsize', 12)
xlabel('$t$ (days)' , 'Interpreter', 'latex', 'FontSize', 12)

subplot(3,1,3)
hold on
plot(0:tBreaks(end), Ddata(1:tBreaks(end)+1), 'o', 'LineWidth', 2, 'Color', colorlist(3,:), 'MarkerSize', 5);
plot(0:tBreaks(end), D_IM(1:tBreaks(end)+1), 'LineWidth', 2, 'Color', colorlist(2,:));
plot(0:tBreaks(end), D_IGM(1:tBreaks(end)+1), 'LineWidth', 2, 'Color', colorlist(1,:));

xlim([0 539])
legend('Real data', 'Sim. IM', 'Sim. IGM', 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'northwest')
ylabel('$D(t)$', 'Interpreter', 'latex', 'FontSize', 12)
set(gca, 'TickLabelInterpreter', 'latex', 'fontsize', 12)
xlabel('$t$ (days)' , 'Interpreter', 'latex', 'FontSize', 12)

%% Plot parameters

beta0_IGM = [];
gamma_IGM = [];
tau_IGM = [];
omega_IGM = [];
zeta_IGM = [];
sigmaC0_IGM = [];
sigmaN0_IGM = [];
alpha_IGM = [];
e_IGM = [];
a_IGM = [];

beta0_IM = [];
gamma_IM = [];
tau_IM = [];
omega_IM = [];
zeta_IM = [];
sigmaC0_IM = [];
sigmaN0_IM = [];
alpha_IM = [];
e_IM = [];
a_IM = [];

for i=1:W % Extract the estimated parameter values
    
    val = getParameter(par_obj_IGM, W, thetaOptim_IGM, 'beta0', i);
    beta0_IGM = [beta0_IGM; val*ones(Tw, 1)];
    
    val = getParameter(par_obj_IGM, W, thetaOptim_IGM, 'gamma', i);
    gamma_IGM = [gamma_IGM; val*ones(Tw, 1)];
    
    val = getParameter(par_obj_IGM, W, thetaOptim_IGM, 'tau', i);
    tau_IGM = [tau_IGM; val*ones(Tw, 1)];
    
    val = getParameter(par_obj_IGM, W, thetaOptim_IGM, 'omega', i);
    omega_IGM = [omega_IGM; val*ones(Tw, 1)];
    
    val = getParameter(par_obj_IGM, W, thetaOptim_IGM, 'zeta', i);
    zeta_IGM = [zeta_IGM; val*ones(Tw, 1)]; 
    
    val = getParameter(par_obj_IGM, W, thetaOptim_IGM, 'sigmaC0', i);
    sigmaC0_IGM = [sigmaC0_IGM; val*ones(Tw, 1)];   
    
    val = getParameter(par_obj_IGM, W, thetaOptim_IGM, 'sigmaN0', i);
    sigmaN0_IGM = [sigmaN0_IGM; val*ones(Tw, 1)];   
    
    val = getParameter(par_obj_IGM, W, thetaOptim_IGM, 'alpha', i);
    alpha_IGM = [alpha_IGM; val*ones(Tw, 1)];  
    
    val = getParameter(par_obj_IGM, W, thetaOptim_IGM, 'e', i);
    e_IGM = [e_IGM; val*ones(Tw, 1)];     
    
    val = getParameter(par_obj_IGM, W, thetaOptim_IGM, 'a', i);
    a_IGM = [a_IGM; val*ones(Tw, 1)]; 
    
    %% IM
    val = getParameter(par_obj_IM, W, thetaOptim_IM, 'beta0', i);
    beta0_IM = [beta0_IM; val*ones(Tw, 1)];
    
    val = getParameter(par_obj_IM, W, thetaOptim_IM, 'gamma', i);
    gamma_IM = [gamma_IM; val*ones(Tw, 1)];
    
    val = getParameter(par_obj_IM, W, thetaOptim_IM, 'tau', i);
    tau_IM = [tau_IM; val*ones(Tw, 1)];
    
    val = getParameter(par_obj_IM, W, thetaOptim_IM, 'omega', i);
    omega_IM = [omega_IM; val*ones(Tw, 1)];
    
    val = getParameter(par_obj_IM, W, thetaOptim_IM, 'zeta', i);
    zeta_IM = [zeta_IM; val*ones(Tw, 1)]; 
    
    val = getParameter(par_obj_IM, W, thetaOptim_IM, 'sigmaC0', i);
    sigmaC0_IM = [sigmaC0_IM; val*ones(Tw, 1)];   
    
    val = getParameter(par_obj_IM, W, thetaOptim_IM, 'sigmaN0', i);
    sigmaN0_IM = [sigmaN0_IM; val*ones(Tw, 1)];   
    
    val = getParameter(par_obj_IM, W, thetaOptim_IM, 'alpha', i);
    alpha_IM = [alpha_IM; val*ones(Tw, 1)];  
    
    val = getParameter(par_obj_IM, W, thetaOptim_IM, 'e', i);
    e_IM = [e_IM; val*ones(Tw, 1)];     
    
    val = getParameter(par_obj_IM, W, thetaOptim_IM, 'a', i);
    a_IM = [a_IM; val*ones(Tw, 1)];     
end


%% alpha, omega, zeta

figure
subplot(1,3,1)
hold on
plot(0:tBreaks(end), alpha_IM, 'LineWidth', 2, 'Color', colorlist(2,:));
plot(0:tBreaks(end), alpha_IGM, 'LineWidth', 2, 'Color', colorlist(1,:));

xlim([0 539]) 
legend('IM', 'IGM', 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'northwest')
ylabel('$\alpha(t)$', 'Interpreter', 'latex', 'FontSize', 12)
set(gca, 'TickLabelInterpreter', 'latex', 'fontsize', 12)
xlabel('$t$ (days)' , 'Interpreter', 'latex', 'FontSize', 12)

subplot(1,3,2)
hold on
plot(0:tBreaks(end), omega_IM, 'LineWidth', 2, 'Color', colorlist(2,:));
plot(0:tBreaks(end), omega_IGM, 'LineWidth', 2, 'Color', colorlist(1,:));

xlim([0 539])
legend('IM', 'IGM', 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'northwest')
ylabel('$\omega(t)$', 'Interpreter', 'latex', 'FontSize', 12)
set(gca, 'TickLabelInterpreter', 'latex', 'fontsize', 12)
xlabel('$t$ (days)' , 'Interpreter', 'latex', 'FontSize', 12)

subplot(1,3,3)
hold on
plot(0:tBreaks(end), zeta_IM, 'LineWidth', 2, 'Color', colorlist(2,:));
plot(0:tBreaks(end), zeta_IGM, 'LineWidth', 2, 'Color', colorlist(1,:));

xlim([0 539]) 
legend('IM', 'IGM', 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'northwest')
ylabel('$\zeta(t)$', 'Interpreter', 'latex', 'FontSize', 12)
set(gca, 'TickLabelInterpreter', 'latex', 'fontsize', 12)
xlabel('$t$ (days)' , 'Interpreter', 'latex', 'FontSize', 12)

%% a, G(I,D) - e, x, ex - beta

figure
subplot(1,3,1)
hold on
yyaxis left
plot(0:tBreaks(end), a_IGM, 'LineWidth', 2, 'Color', colorlist(1,:));
set(gca, 'YColor', colorlist(1,:));

yyaxis right
G_IGM = (I_IGM+D_IGM-a_IGM)/M;
plot(0:tBreaks(end), G_IGM, 'LineWidth', 2, 'Color', colorlist(4,:));
plot(0:tBreaks(end), 0*G_IGM, '--', 'LineWidth', 2, 'Color',  [0 0 0]);
set(gca, 'YColor', colorlist(4,:));

xlim([0 539])
legend('$a$', '$G(I,D)$', 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'northwest')
set(gca, 'TickLabelInterpreter', 'latex', 'fontsize', 12)
xlabel('$t$ (days)' , 'Interpreter', 'latex', 'FontSize', 12)

subplot(1,3,2)
hold on
plot(0:tBreaks(end), e_IGM, 'LineWidth', 2, 'Color', colorlist(1,:));
plot(0:tBreaks(end), x_IGM, 'LineWidth', 2, 'Color', colorlist(4,:));
plot(0:tBreaks(end), e_IGM.*x_IGM, 'LineWidth', 2, 'Color', colorlist(5,:));

xlim([0 539]) 
legend('$e$', '$x$', '$ex$', 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'northwest')
set(gca, 'TickLabelInterpreter', 'latex', 'fontsize', 12)
xlabel('$t$ (days)' , 'Interpreter', 'latex', 'FontSize', 12)

subplot(1,3,3)
hold on
plot(0:tBreaks(end), beta0_IM, 'LineWidth', 2, 'Color', colorlist(2,:));
plot(0:tBreaks(end), beta0_IGM.*(1-e_IGM.*x_IGM), 'LineWidth', 2, 'Color', colorlist(1,:));

xlim([0 539]) 
legend('$\beta$ (IM)', '$\beta$ (IGM)', 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'northwest')
set(gca, 'TickLabelInterpreter', 'latex', 'fontsize', 12)
xlabel('$t$ (days)' , 'Interpreter', 'latex', 'FontSize', 12)

return

%% U
figure
% subplot(1,2,1)
hold on
yyaxis left
% plot(0:tBreaks(end), res_nox.Y(1:tBreaks(end)+1), 'LineWidth', 2);
plot(0:tBreaks(end), U_IGM(1:tBreaks(end)+1), 'LineWidth', 2, 'Color', colorlist(1,:));
ylabel('$U(t)$', 'Interpreter', 'latex', 'FontSize', 12)
set(gca, 'YColor', colorlist(1,:));

yyaxis right


data_yex = 100*U_IGM(1:tBreaks(end)+1)./(U_IGM(1:tBreaks(end)+1)+I_IGM(1:tBreaks(end)+1));
% data_yex_2 = 100*mypos(res_yex.Y(1:tBreaks(end)+1))./mypos(res_yex.C(1:tBreaks(end)+1));
plot(0:tBreaks(end), data_yex, 'LineWidth', 2, 'Color', colorlist(4,:));
% plot(0:tBreaks(end), data_yex_2, 'LineWidth', 2);       
ylabel('$U(t)$ (\%)', 'Interpreter', 'latex', 'FontSize', 12)
set(gca, 'YColor', colorlist(4,:));
set(gca, 'TickLabelInterpreter', 'latex', 'fontsize', 12)
xlim([0 539])
xlabel('$t$ (days)' , 'Interpreter', 'latex', 'FontSize', 12)


%% R0

RO_IM = beta0_IM./(gamma_IM+tau_IM);                                                   
RO_IGM = beta0_IGM.*(1-e_IGM.*x_IGM)./(gamma_IGM+tau_IGM);
RO_IGM_0 = beta0_IGM./(gamma_IGM+tau_IGM);
RO_IGM_ott = beta0_IGM.*(1-e_IGM)./(gamma_IGM+tau_IGM);
R0_diff_data = zeros(numel(Idata)-1, 1);
for i=1:numel(Idata)-1
    R0_diff_data(i) = 1 + ((Idata(i+1)-Idata(i))/Idata(i))/(1/14 + 0.166);
end

k1p_est = M*(1-1/RO_IGM_0(1));
fprintf('k1 primo: %.4f\n', k1p_est);


% R0_diff_data(448) = (R0_diff_data(449) + R0_diff_data(447))/2;

figure; %('Position', [285         401        1270         270])
% subplot(1,3,1)
hold on

x_patch = 0:tBreaks(end);
x_patch = [0:tBreaks(end), x_patch(end:-1:1)];

y_patch_1 = [RO_IGM_0; RO_IGM(end:-1:1)];
hp1 = patch(x_patch, y_patch_1, 'c');
hp1.FaceColor = colorlist(4,:);
hp1.FaceAlpha = 0.3;
y_patch_2 = [RO_IGM; RO_IGM_ott(end:-1:1)];
hp2 = patch(x_patch, y_patch_2, 'g');
hp2.FaceColor = colorlist(5,:);
hp2.FaceAlpha = 0.3;

h = zeros(4, 1);
h(1) = plot(1:tBreaks(end), R0_diff_data, 'o', 'LineWidth', 2, 'Color', colorlist(3,:), 'MarkerSize', 3);

h(2) = plot(0:tBreaks(end), RO_IGM, 'LineWidth', 2, 'Color', colorlist(1,:));

h(3) = plot(0:tBreaks(end), RO_IGM_0, 'LineWidth', 2, 'Color', colorlist(4,:));

h(4) = plot(0:tBreaks(end), RO_IGM_ott, 'LineWidth', 2, 'Color', colorlist(5,:));


plot(0:tBreaks(end), 1+0*RO_IGM_ott, '--', 'LineWidth', 2, 'Color', [0 0 0 ]);

% 
% yyaxis right
% h(5) = plot(2:tBreaks(end), diff(R0_diff_data), 'r', 'LineWidth', 2);
% set(gca, 'YColor', 'r')

% plot(0:tBreaks(end), alpha_IGM, 'LineWidth', 2, 'Color', colorlist(1,:));

xlim([0 539]) %tBreaks(end)])

legend(h, 'Real data', '$\mathcal{R}(t)$ (IGM)', '$\mathcal{R}_0$ (IGM)', '$\underline{\mathcal{R}}$ (IGM)', 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'northwest')
% ylabel('$\alpha(t)$', 'Interpreter', 'latex', 'FontSize', 12)
set(gca, 'TickLabelInterpreter', 'latex', 'fontsize', 12)
xlabel('$t$ (days)' , 'Interpreter', 'latex', 'FontSize', 12)

yl = ylim;
ylim([0.2 1.5])



%%
figure
hs = scatter(e_IGM(1:30:540), I_IGM(1:30:540), 200, 1:18, 'filled');
hold on
hs2 = scatter(e_IGM(1:30:510), I_IGM(31:30:540), 200, 1:17, 'filled');
hs2.Marker = 'd';
for i=1:18
    text(e_IGM(1+(i-1)*30), I_IGM(1+(i-1)*30), num2str(i));
end



% AA = zeros(17);
% for i=1:16
%     AA(i,i+1) = 1;
% end
% G = digraph(AA);


figure
% subplot(1,2,1)

myranges_I = {7, [5 6], [1 2 3 4 8:17]};
myranges_D = {[6 7], [], [1:5  8:17]};

for j=1:3

    if (j < 3)
        subplot(6,1,j)
    else
        subplot(6,1,[3 4 5 6])
    end
        
    i_tmp = zeros(18, 1);
    d_tmp = zeros(18, 1);
%     i_diff = diff(I_IGM(1:540));
    d_diff = diff(D_IGM(1:540));

    for i=1:18
        
        if (i<18)
            d_tmp(i) = max(d_diff(1+(i-1)*30:30*i));
            i_tmp(i) = max(I_IGM(1+(i-1)*30:30*i));
        else
            d_tmp(i) = max(d_diff(1+(i-1)*30:30*i-1));
            i_tmp(i) = max(I_IGM(1+(i-1)*30:30*i-1));
        end
    end
      

    y =100*(i_tmp(2:18)-i_tmp(1:17))./i_tmp(1:17);
    y(12) = -5;
    
    
    x = e_IGM(1:30:510);
    hold on
    plot([0, 1], 0*[0 1], 'k--')
    hs = scatter(x,y,300,1:17,'filled');
    hs.MarkerFaceColor = colorlist(3,:);
    hs.MarkerFaceAlpha = 0.7;



    for i=myranges_I{j} %1:17
        htxt = text(x(i)+0.02, y(i)+0.02, num2str(i));
        htxt.FontSize = 12;
    end
%     ylabel('$I(t)$', 'Interpreter', 'latex', 'FontSize', 12)
    set(gca, 'TickLabelInterpreter', 'latex', 'fontsize', 16)
    

    y =100*(d_tmp(2:18)-d_tmp(1:17))./d_tmp(1:17);
    y(12) = -5;
    
    x = e_IGM(1:30:510);
    hold on
%     plot([min(x), max(x)], 0*[min(x), max(x)], 'k--')


    hs2 = scatter(x,y,300,1:17,'filled');
    hs2.MarkerFaceColor = colorlist(4,:);
    hs2.MarkerFaceAlpha = 0.7;
    for i=myranges_D{j} %1:17
        htxt = text(x(i)-0.02, y(i)-0.02, num2str(i));
        htxt.FontSize = 12;
    end
%     ylabel('$\dot{D}(t)$', 'Interpreter', 'latex', 'FontSize', 12)
    set(gca, 'TickLabelInterpreter', 'latex', 'fontsize', 16)
    
    
    if (j == 1)
        ylim([450 560])

        set(gca,'XTick', [], 'YTick', 500)        
    elseif (j==2)
        ylim([180 205])
        set(gca,'XTick', [], 'YTick', 195)
    else
        ylim([-95 35]);
        xlabel('$e$' , 'Interpreter', 'latex', 'FontSize', 16)
        ylabel('$\Delta \%$' , 'Interpreter', 'latex', 'FontSize', 16)
        legend([hs hs2], '$I$', '$D$',  'Interpreter', 'latex', 'FontSize', 16)
    end    
    
    xlim([0 1])

end    
    
return

% % %  REGIONI - TRIMESTRI


I_IGM_tmp_end = zeros(6, 1);
E_IGM_tmp_end = zeros(6, 1);
A_IGM_tmp_end = zeros(6, 1);
for i=1:6
    I_IGM_tmp_end(i) = max(I_IGM(1+(i-1)*90:i*90));
    E_IGM_tmp_end(i) = mean(e_IGM(1+(i-1)*90:i*90));
    A_IGM_tmp_end(i) = mean(a_IGM(1+(i-1)*90:i*90));    
end

figure
hold on
scatter(E_IGM_tmp_end, A_IGM_tmp_end, 200, I_IGM_tmp_end, 'filled')

for i=1:6
    ht = text(E_IGM_tmp_end(i), A_IGM_tmp_end(i), sprintf('%d', i));
    set(ht, 'FontSize', 16)
end
aut = autumn;
aut = aut(end:-1:1,:);
colormap(aut)
colorbar

% % % REGIONI - PRED UN MESE AVANTI PER DIVERSI e ed a

e_grid = 0.1:0.1:0.9;
a_grid = linspace(1e5, 4e5, 9);

e_pred = zeros(numel(e_grid)*numel(a_grid), 1);
a_pred = zeros(numel(e_grid)*numel(a_grid), 1);
I_pred = zeros(numel(e_grid)*numel(a_grid), 1);

% Idata, ...
% Rdata, ...
% Ddata, ...

tBreaks = zeros(W-1, 1);
tBreaks(1) = Tw-1;
for i=2:W-1
    tBreaks(i) = tBreaks(i-1) + Tw;
end

h = 0;
for i=1:numel(e_grid)
    for j=1:numel(a_grid)

        [i j]
        
        
        thetaOptim_pred_IGM = [thetaOptim_IGM; thetaOptim_IGM(end-(n_tvar_IM-1):end)]; %thetaOptim_IGM(end-(n_tvar_IM-1):end)];
        

        thetaOptim_pred_IGM(7+4+(19-1)*5) = e_grid(i);
        thetaOptim_pred_IGM(7+5+(19-1)*5) = a_grid(j);  
%         
%         thetaOptim_pred_IGM(7+4+(20-1)*5) = e_grid(i);
%         thetaOptim_pred_IGM(7+5+(20-1)*5) = a_grid(j);          
    
        [S_IGM_tmp, U_IGM_tmp, I_IGM_tmp, R_IGM_tmp, D_IGM_tmp, x_IGM_tmp] = simode(thetaOptim_pred_IGM, ...
                                                                        Idata(1), ...
                                                                        Rdata(1), ...
                                                                        Ddata(1), ...
                                                                        M, ...
                                                                        tBreaks, ...
                                                                        par_obj_IGM);        
                                                                    
        h = h + 1;
        I_pred(h) = I_IGM_tmp(end);
        e_pred(h) = e_grid(i);
        a_pred(h) = a_grid(j);
        
    end
end

scatter(e_pred, a_pred, 100, I_pred, 'Marker', 'd') 


figure
scatter(e_pred, a_pred, 100, I_pred, 'filled') 

aut = autumn;
aut = aut(end:-1:1,:);
colormap(aut)
colorbar






