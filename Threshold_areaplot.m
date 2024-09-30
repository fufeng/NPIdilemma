% Parameters
beta = 1; % Example fixed value for beta
C_FC = 1; % Example fixed value for C_FC
C_I = 10; % Example fixed value for C_I
C_SD_range = linspace(C_FC, C_I, 1000); % Define a range for C_SD, ensuring it's within valid bounds

% Calculate beta_FC* for the range of C_SD
Beta_FC_star = beta * (log(1 - (C_SD_range - C_FC) / C_I) ./ log(1 - C_SD_range / C_I));


% Define the RGB values for the colors
colors = [
    255, 192, 203;  % Light Pink
    255, 182, 193;  % Light Pink 1
    255, 105, 180;  % Hot Pink
    255, 20, 147;   % Deep Pink
    219, 112, 147;  % Pale Violet Red
    186, 85, 211;   % Medium Orchid
    147, 112, 219;  % Medium Purple
    138, 43, 226;   % Blue Violet
    139, 0, 139;    % Dark Magenta
    72, 61, 139     % Dark Slate Blue
] / 255; % Normalize RGB values to [0, 1]

rose = [0.961, 0.51,0.933]; 
lightBlue = [0.298, 0.69, 0.961]; 
purple = [1, 1,0.5]; 
pastelGreen = [0.298, 0.85, 0.851]; 

% Plotting
figure;
hold on;

C_SD_range = C_SD_range/C_I;
% Fill area above and below the curve in different colors
fill([C_SD_range fliplr(C_SD_range)], [Beta_FC_star fliplr(1.2*ones(1,length(C_SD_range)))], pastelGreen, 'LineStyle', 'none','DisplayName','(iii) I_2* < I_1*, no I_3* or I_4*');
fill([C_SD_range(1:646) fliplr(C_SD_range(1:646))], [Beta_FC_star(1:646) fliplr(0.7613*ones(1,length(C_SD_range(1:646))))], lightBlue, 'LineStyle', 'none','DisplayName','(i) I_2* < I_1* < I_3* < I_4*');
fill([C_SD_range(646:end) fliplr(C_SD_range(646:end))], [Beta_FC_star(646:end) fliplr(0.7613*ones(1,length(C_SD_range(646:end))))], purple, 'LineStyle', 'none','DisplayName','(iv) I_3* < I_4* < I_2* < I_1*');
fill([C_SD_range fliplr(C_SD_range)], [Beta_FC_star fliplr(0*ones(1,length(C_SD_range)))], rose, 'LineStyle', 'none','DisplayName','(ii) I_3* < I_2* < I_1* < I_4*'); % Above

xline(6.8108/C_I,'b--','LineWidth',5,'DisplayName','Scenario separation')
%plot(C_SD_range, Beta_FC_star, 'b--', 'LineWidth', 5, 'DisplayName',''); % Plot beta_FC* as a function of C_SD


xlabel('Relative cost of social distancing, C_{SD}/C_I');
ylabel('Relative risk of infection (FC), \beta_{FC}/\beta');
%title('Thresholds for different behaviors');
%ylim(yLimits); % Reset y-axis limits in case filling altered them
ylim([0 1]);
%xlim([1 10]);
xlim([0.1 1]);
set(gca,'FontSize',20)
legend
hold off;
view([90 -90]);