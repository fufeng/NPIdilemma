% Parameters
C_FC = 2;
C_I = 5;
beta = 10;
beta_FC_values = [0.2*beta, 0.315*beta, 0.4*beta]; % Different values of beta_FC to illustrate different scenarios
line_style = ['- ';'-.';'--'];

% Range of I3* values for plotting
I3_star_range = linspace(0, 1, 1000);

% Create a figure
figure;
hold on;

% Plot function for different beta_FC values
for i = 1:length(beta_FC_values)
    beta_FC = beta_FC_values(i);
    func_values = exp(-beta_FC * I3_star_range) - exp(-beta * I3_star_range) - C_FC/C_I;
    plot(I3_star_range, func_values,line_style(i,:), 'DisplayName', ['\beta_{FC} = ', num2str(beta_FC)],'LineWidth',5);
end

mycolors = [0.655 0.824 0.91; 0.961 0.545 0.871; 0.424 0.388 0.8];
ax = gca; 
ax.ColorOrder = mycolors;

% Add plot formatting
xlabel('Disease prevalence, I');
ylabel('Payoff difference, \pi_{FC} - \pi_{NSD}');
xline(0.081,'-.','Linewidth',3,'DisplayName','I_3*')
xline(0.437,':','Linewidth',3,'DisplayName','I_4*')
%title('Function of I_3^* for Different \beta_{FC} Values');
set(gca,'fontsize',20)
legend('show');
yline(0,'--','LineWidth',3, 'HandleVisibility', 'off')
grid off;
hold off;
