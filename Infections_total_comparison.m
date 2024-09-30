% Graphics for different values CI
% Set parameters
beta = 10; beta_fc = 4; gamma = 4; omega = 3; kappa = 3;
C_SD = 2; C_I = 10; C_FC = 1;

tspan = [0, 1000];

Y0 = [0.99, 0.01, 0, 0.8];

rose = [0.961, 0.51,0.933]; 
lightBlue = [0.298, 0.69, 0.961]; 
purple = [1, 1,0.5]; 
pastelGreen = [0.298, 0.85, 0.851]; 

% Solve ODE for different CI
range= 0:0.01*beta:beta;
r_infty1 = ones(1,length(range));
r_infty2 = ones(1,length(range));
r_infty3 = ones(1,length(range));
I_1_star = ones(1,length(range));
perf_ad1 = ones(1,length(range));
perf_ad2 = ones(1,length(range));

i=0;
for beta_fc = range
    %I_star = -log(1 - C_sd/ci)/beta;
    %[~,SIR_SD1] = ode45(@(t,Y) SIR_SD_FC(t, Y, beta, beta_fc, gamma, omega, kappa, 2, C_I, C_FC), tspan, Y0);
    %[~,SIR_SD2] = ode45(@(t,Y) SIR_SD_FC(t, Y, beta, beta_fc, gamma, omega, kappa, 3, C_I, C_FC), tspan, Y0);
    [~,SIR_SD3] = ode23(@(t,Y) SIR_SD_FC(t, Y, beta, beta_fc, gamma, omega, kappa, 15, C_I, C_FC), tspan, Y0);
    i = i+1;
    %r_infty1(i) = 1-SIR_SD1(end,1)-SIR_SD1(end,2);
    %r_infty2(i) = 1-SIR_SD2(end,1)-SIR_SD2(end,2);
    r_infty3(i) = 1-SIR_SD3(end,1)-SIR_SD3(end,2);

    I_1_star1 = 10;%-1/beta*log(1-C_SD/C_I);
    I_1_star2 = 10; %-1/beta*log(1-10/C_I);

    func = @(I3_star) exp(-beta_fc * I3_star) - exp(-beta * I3_star) - C_FC/C_I;

    roots_1 = 10;
    
    % Attempt to find a root with initial guess near 0
    options = optimset('Display', 'off'); % Turn off fzero output
    [root1, fval1, exitflag1] = fzero(func, 0, options);
    if exitflag1 > 0 % Check if fzero succeeded
        roots_1 = root1;
    end

    if roots_1 < 10
        I_3_star1 = roots_1;
    else
        I_3_star1 = 10;
    end
    display(I_3_star1)

    func = @(I3_star) exp(-beta_fc * I3_star) - exp(-beta * I3_star) - 15/C_I;

    roots_1 = 10;

        % Attempt to find a root with initial guess near 0.2
    options = optimset('Display', 'off'); % Turn off fzero output
    [root1, fval1, exitflag1] = fzero(func, 0, options);
    if exitflag1 > 0 % Check if fzero succeeded
        roots_1 = root1;
    end

    if roots_1 < 10
        I_3_star2 = roots_1;
    else
        I_3_star2 = 10;
    end



    if I_3_star1 < I_1_star1
        perf_ad1(i) = gamma/beta*lambertw(-exp(-I_3_star1*beta/gamma-1))+1;
    else
        perf_ad1(i) = gamma/beta*lambertw(-exp(-I_1_star1*beta/gamma-1))+1;
    end

    if I_3_star2 < I_1_star1
        perf_ad2(i) = gamma/beta*lambertw(-exp(-I_3_star2*beta/gamma-1))+1;
    else
        perf_ad2(i) = gamma/beta*lambertw(-exp(-I_1_star1*beta/gamma-1))+1;
    end


end

Y0 = [0.99, 0.01, 0.0, 1];

[~,SIR] = ode45(@(t,Y) SIR_SD_FC(t, Y, beta, beta_fc, gamma, omega, kappa, 2, C_I, C_FC), tspan, Y0);

hold on;
% Plot
plot(range(1:77)/beta,perf_ad1(1:77),':',range/beta,r_infty3,'-','Linewidth',4)
yline(1-SIR(end,1)-SIR(end,2),'--k','LineWidth',4)
plot(range(78:end)/beta,ones(length(range(78:end)))*(1-SIR(end,1)-SIR(end,2)),':','Linewidth',4,'color',[0.961 0.545 0.871])

xlabel('Relative effectiveness of face covering \beta_{FC}','FontSize',15); ylabel('Total Infections R(\infty)','FontSize',15);
%axis([10 450 0.35 0.7]);
legend('Perf. Ad. (only FC)','SIR-FC model','SIR model baseline','FontSize',15)
set(gca,'FontSize',15)
%axis([1 26 0.6 0.9])
mycolors = [0.961 0.545 0.871; 0.455 0.624 0.91; 0.298, 0.961, 0.851];
ax = gca; 
ax.ColorOrder = mycolors;

hold off

