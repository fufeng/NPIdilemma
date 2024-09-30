% Parameters
C_FC = 2;
C_I = 5;
beta = 10;
beta_FC_values = linspace(1.5, beta, 10000); % Different values of beta_FC to illustrate different scenarios
beta_FC_plot = linspace(1.5/beta, 1, 10000);

% Initialize arrays to store roots
roots1 = [];
roots2 = [];
beta_FC_for_roots1 = [];
beta_FC_for_roots2 = [];

% Loop over beta_FC values to find roots
for beta_FC = beta_FC_values
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
end

common_beta_FC = linspace(min([beta_FC_for_roots1, beta_FC_for_roots2]), max([beta_FC_for_roots1, beta_FC_for_roots2]), 1000);
common_beta_FC_plot = linspace(min([beta_FC_for_roots1/beta, beta_FC_for_roots2/beta]), max([beta_FC_for_roots1/beta, beta_FC_for_roots2/beta]), 1000);
% Interpolate roots1 and roots2 to this common_beta_FC
interp_roots1 = interp1(beta_FC_for_roots1, roots1, common_beta_FC, 'linear', 'extrap');
interp_roots2 = interp1(beta_FC_for_roots2, roots2, common_beta_FC, 'linear', 'extrap');

% Plotting
figure;
hold on;

% Define vertical boundaries for the plot, adjust these values as needed
ymin = min([interp_roots1, interp_roots2]) - 0.1; % Extend 0.1 below lowest root
ymax = max([interp_roots1, interp_roots2]) + 0.1; % Extend 0.1 above highest root

% Fill area between roots
fill_area_x = [common_beta_FC_plot, fliplr(common_beta_FC_plot)];
fill_area_y = [interp_roots1, fliplr(interp_roots2)];
fill(fill_area_x, fill_area_y, [0.8, 0, 0.4], 'LineStyle', 'none', 'DisplayName', '\pi_{FC} > \pi_{NSD}'); % Fill the area

% Fill the remaining areas in light blue
% Below interp_roots1
fill([common_beta_FC_plot, fliplr(common_beta_FC_plot)], [ymin*ones(1, length(common_beta_FC)), fliplr(interp_roots1)], [0.651, 0.827, 0.988], 'LineStyle', 'none', 'DisplayName', '\pi_{FC} < \pi_{NSD}');
% Above interp_roots2
fill([common_beta_FC_plot, fliplr(common_beta_FC_plot)], [ymax*ones(1, length(common_beta_FC)), fliplr(interp_roots2)], [0.651, 0.827, 0.988], 'LineStyle', 'none', 'HandleVisibility', 'off');

% Remaining area
remaining_area = linspace(common_beta_FC_plot(end),3.6,2);
fill([remaining_area, fliplr(remaining_area)], [ymax*ones(1, length(remaining_area)), fliplr(ymin*ones(1, length(remaining_area)))], [0.651, 0.827, 0.988], 'LineStyle', 'none', 'HandleVisibility', 'off')

% Plot roots1 and roots2 lines for visibility
plot(common_beta_FC_plot, interp_roots1, 'b--', 'DisplayName', '\pi_{FC} = \pi_{NSD}', 'LineWidth', 5);
plot(common_beta_FC_plot, interp_roots2, 'b--', 'LineWidth', 5, 'HandleVisibility', 'off');

xlabel('Relative risk of infection (FC), \beta_{FC}/\beta');
ylabel('Disease prevalence, I');
xlim([1.5/beta 3.5/beta])
ylim([0 0.7])
%title('Payoff of NSD vs FC');
legend('show');
set(gca, 'FontSize', 20);
grid on; % Or grid off, based on your preference
hold off;
