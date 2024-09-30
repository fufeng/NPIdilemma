% Parameters
beta = 10; beta_fc = 4; gamma = 4; omega = 30; kappa = 30;
C_SD = 2; C_I = 10; C_FC = 1;

% Initial conditions [S, I, R, SD]
Y0 = [0.99, 0.01, 0.1, 0.8]; % Initial population distribution

% Time span
tspan = [0 20]; % From day 0 to day 100

% Solve the ODE
[t, Y] = ode23(@(t,Y) SIR_SD_FC(t, Y, beta, beta_fc, gamma, omega, kappa, C_SD, C_I, C_FC), tspan, Y0);

% Colors
rose = [0.961, 0.51, 0.933];
lightBlue = [0.298, 0.69, 0.961];
purple = [0.58, 0.298, 0.961];
pastelGreen = [0.298, 0.961, 0.851];

% Initialize root storage
roots_1 = [];
roots_2 = [];
beta_FC_for_roots1 = [];
beta_FC_for_roots2 = [];

func = @(I3_star) exp(-beta_fc * I3_star) - exp(-beta * I3_star) - C_FC/C_I;

roots_1 = 0;

% Attempt to find a root with initial guess near 0.2
options = optimset('Display', 'off'); % Turn off fzero output
[root1, fval1, exitflag1] = fzero(func, 0, options);
if exitflag1 > 0 % Check if fzero succeeded
    roots_1 = root1;
    beta_FC_for_roots1(end+1) = beta_fc;
end

[root2, fval2, exitflag2] = fzero(func, 0.5, options);
if exitflag2 > 0 && abs(root2 - root1) > 1e-3 % Ensure it's a different root
    roots_2 = root2;
    beta_FC_for_roots2(end+1) = beta_fc;
end

% Plot the results
figure;
hold on;
%plot(t, Y(:,1), 'Color', rose, 'LineWidth', 4);
plot(t, Y(:,2), 'Color', lightBlue, 'LineStyle', '-','LineWidth', 4);
%plot(t, Y(:,3), 'Color', purple, 'LineStyle', ':','LineWidth', 4);
%plot(t, (1-Y(:,3)-Y(:,4)), 'Color', pastelGreen, 'LineStyle', '--', 'LineWidth', 4);
yline(-1/beta_fc*log(1-(C_SD-C_FC)/C_I),'-.k','LineWidth',3)
yline(-1/beta*log(1-C_SD/C_I),'--r','LineWidth',3)
if roots_1 > 0
    yline(roots_1,'-.b','LineWidth',3)
end

hold off;

legend('Infected','I_2^*', 'I_1^*','I_3^*');
xlabel('Time');
ylabel('Population Fraction');
%title('SIR Model with Social Distancing and Face Covering');
set(gca,'Fontsize',20)


figure;
plot(t, Y(:,3), 'Color', purple, 'LineStyle', ':','LineWidth', 4);
hold on;
plot(t, (1-Y(:,3)-Y(:,4)), 'Color', pastelGreen, 'LineStyle', '--', 'LineWidth', 4);
legend('SD','FC');
xlabel('Time');
ylabel('Population Fraction');
%title('SIR Model with Social Distancing and Face Covering');
set(gca,'Fontsize',20)



vertex_labels = {'SD', 'FC', 'NSD'};
% Normalize each set of points such that Y1 + Y2 + Y3 = 1
Y_norm = [Y(:,3),1- Y(:,3)- Y(:,4), Y(:,4)];

% Convert the normalized data to barycentric coordinates for plotting on the simplex
x_simplex = 0.5 * (2 * Y_norm(:,2) + Y_norm(:,3)) ./ (Y_norm(:,1) + Y_norm(:,2) + Y_norm(:,3));
y_simplex = (sqrt(3)/2) * Y_norm(:,3) ./ (Y_norm(:,1) + Y_norm(:,2) + Y_norm(:,3));

% Create a simplex (triangle) for the plot
figure;
hold on;
plot([0 1 0.5 0], [0 0 sqrt(3)/2 0], 'k-', 'LineWidth', 2); % The triangle boundary

% Plot the solutions over time in the simplex
scatter(x_simplex, y_simplex, 50, 1:length(x_simplex), 'filled'); % Color by time

quiver(x_simplex(1:end-1), y_simplex(1:end-1), ...
       x_simplex(2:end) - x_simplex(1:end-1), ...
       y_simplex(2:end) - y_simplex(1:end-1), ...
       0, 'k','LineWidth', 2, 'MaxHeadSize', 2); % 0 scaling means actual size, 'k' is black color for arrows

% Label the vertices of the triangle
text(0, 0, vertex_labels{1}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
text(1, 0, vertex_labels{2}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
text(0.5, sqrt(3)/2, vertex_labels{3}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');


% Colorbar to indicate time progression
hColorbar = colorbar;
% Customize tick labels
tickLabels = {'Start', 'End'};
set(hColorbar, 'YTick', [min(hColorbar.Limits) max(hColorbar.Limits)], 'YTickLabel', tickLabels);

% Add a label at the top of the colorbar
% ylabel(hColorbar, 'Time', 'Rotation', 0, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontWeight', 'bold');

axis equal;
axis off;
hold off;
saveas(gcf, 'simplexplot.svg');




% Define video writer object
video = VideoWriter('simplex_animation.mp4', 'MPEG-4'); % Create a .mp4 file
video.FrameRate = 5; % Set the frame rate (adjust as needed)
open(video); % Open the video file for writing

% Create a simplex (triangle) for the plot
figure;
hold on;
xlim([0,1])
ylim([0,1])
plot([0 1 0.5 0], [0 0 sqrt(3)/2 0], 'k-', 'LineWidth', 2); % The triangle boundary
axis equal;

% Label the vertices of the triangle
text(0, 0, vertex_labels{1}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
text(1, 0, vertex_labels{2}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
text(0.5, sqrt(3)/2, vertex_labels{3}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

axis off;

% Initialize plot objects for real-time updating
h_scatter = scatter(x_simplex(1), y_simplex(1), 50, 1, 'filled');
h_quiver = quiver(0, 0, 0, 0, 0, 'k');  % Initially zero-length arrow

% Loop to update the plot in real-time
for i = 2:length(x_simplex)
    % Update scatter point
    set(h_scatter, 'XData', x_simplex(1:i), 'YData', y_simplex(1:i), 'CData', (1:i)/ length(x_simplex));
    
    % Update arrow to indicate direction of time
    set(h_quiver, 'XData', x_simplex(i-1), ...
                  'YData', y_simplex(i-1), ...
                  'UData', x_simplex(i) - x_simplex(i-1), ...
                  'VData', y_simplex(i) - y_simplex(i-1));

    % Capture the frame and write to the video file
    frame = getframe(gcf);
    writeVideo(video, frame);

    % Pause for a short duration to create an animation effect
    pause(0.1);  % Adjust the pause time to speed up or slow down the animation
    
    % Draw the updates
    drawnow;
end

% Finalize the plot
%colorbar;
% xlabel('Barycentric Coordinate X');
% ylabel('Barycentric Coordinate Y');
% title('Real-Time Animation of ODE Solution in Simplex');
axis equal;

hold off;

close(video);
