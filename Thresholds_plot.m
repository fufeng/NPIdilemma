% Parameters
C_FC = 10;
C_SD = 45;
%C_SD = 20;
C_I = 50;
beta = 15;
beta_FC_values = linspace(3, beta, 10000); % Different values of beta_FC to illustrate different scenarios

% Initialize arrays to store roots
roots1 = [];
roots2 = [];
beta_FC_for_roots1 = [];
beta_FC_for_roots2 = [];
I_2 = [];

% Loop over beta_FC values to find roots
for beta_FC = beta_FC_values
    % Roots for I_3/I_4
    func = @(I3_star) exp(-beta_FC * I3_star) - exp(-beta * I3_star) - C_FC/C_I;
    
    % Attempt to find a root with initial guess near 0.2
    options = optimset('Display', 'off'); % Turn off fzero output
    [root1, fval1, exitflag1] = fzero(func, 0, options);
    if exitflag1 > 0 % Check if fzero succeeded
        roots1(end+1) = root1;
        beta_FC_for_roots1(end+1) = beta_FC;
    end
    
    % Attempt to find a different root with initial guess near 0.8
    [root2, fval2, exitflag2] = fzero(func, 0.5, options);
    if exitflag2 > 0 && abs(root2 - root1) > 1e-3 % Ensure it's a different root
        roots2(end+1) = root2;
        beta_FC_for_roots2(end+1) = beta_FC;
    end

    % Root for I_2
    I_2(end+1) = -1/beta_FC*log(1-(C_SD-C_FC)/C_I);
end

I_1 = -1/beta*log(1-C_SD/C_I);

% Find the common values of beta_FC
common_beta_FC = linspace(min([beta_FC_for_roots1, beta_FC_for_roots2]), max([beta_FC_for_roots1, beta_FC_for_roots2]), 1000);

% Interpolate roots1 and roots2 to this common_beta_FC
interp_roots1 = interp1(beta_FC_for_roots1, roots1, common_beta_FC, 'linear', 'extrap');
interp_roots2 = interp1(beta_FC_for_roots2, roots2, common_beta_FC, 'linear', 'extrap');


% Plotting
figure;
hold on;

% Adjusted color scheme for more rose and light blue tones
rose = [0.961, 0.51,0.933]; 
lightBlue = [0.298, 0.69, 0.961]; 
purple = [0.58, 0.298,0.961]; 
pastelGreen = [0.298, 0.961, 0.851]; 

% Assuming these variables are defined as per your scenario
% common_beta_FC, interp_roots1, interp_roots2, beta_FC_values, I_2, I_1

% Plot the thresholds
yline(I_1, ':', 'Color', pastelGreen, 'LineWidth', 4, 'DisplayName', 'Threshold (I_1^*)');
plot(beta_FC_values/beta, I_2, '-.', 'Color', lightBlue, 'LineWidth', 4, 'DisplayName', 'Threshold (I_2^*)');
plot(common_beta_FC/beta, interp_roots1, '-', 'Color', rose, 'DisplayName', 'Threshold (I_3^*)', 'LineWidth', 4);
plot(common_beta_FC/beta, interp_roots2, '--', 'Color', purple, 'DisplayName', 'Threshold (I_4^*)', 'LineWidth', 4);

xlabel('Relative risk of infection (FC), \beta_{FC}/\beta');
ylabel('Disease prevalence, I');
xlim([4/beta 10/beta]);
ylim([0 0.5]);
%title('Payoff Thresholds');
legend('show');
set(gca, 'FontSize', 20);
grid off;
hold off;
