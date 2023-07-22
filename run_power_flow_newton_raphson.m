% Clear the command window
clc;

% Define the variables
x1 = 0;
x2 = 1;
u = 1;

% Define the functions
h1 = @(x1, x2, u) 4 * u * x2 * sin(x1);
h2 = @(x1, x2, u) 4 * x2^2 - 4 * u * x2 * cos(x1);

% Set the constants
b1 = -0.60;
b2 = -0.30;

% Set the convergence threshold
tolerance = 1e-6;

% Iteration counter
iteration = 1;

while true
    % Calculate the differences (mismatches)
    delta_G1 = b1 - h1(x1, x2, u);
    delta_G2 = b2 - h2(x1, x2, u);

    % Set derivatives of G1 & G2 with respect to x1 & x2
    derivative_G1_dx1 = 4 * u * x2 * cos(x1);
    derivative_G1_dx2 = 4 * u * sin(x1);
    derivative_G2_dx1 = -4 * u * x2 * (-sin(x1));
    derivative_G2_dx2 = 8 * x2 - 4 * u * cos(x1);

    % Create the Jacobian matrix
    Jacobian = [derivative_G1_dx1, derivative_G1_dx2; derivative_G2_dx1, derivative_G2_dx2];

    % Print the mismatches and the Jacobian matrix for this iteration
    disp(['Iteration ', num2str(iteration)]);
    disp('Mismatch for G1: delta_G1 =');
    disp(delta_G1);

    disp('Mismatch for G2: delta_G2 =');
    disp(delta_G2);

    disp('Jacobian =');
    disp(Jacobian);

    % Calculate the updates delta_x1 and delta_x2 using the given formula
    delta_G = [delta_G1; delta_G2];
    Jacobian_inverse = inv(Jacobian);
    updates = Jacobian_inverse * delta_G;

    % Extract the updates for x1 and x2
    delta_x1 = updates(1);
    delta_x2 = updates(2);

    % Update the variables
    x1 = x1 + delta_x1;
    x2 = x2 + delta_x2;

    % Display the updates for this iteration
    disp(['Updates:']);
    disp(['delta_x1 = ', num2str(delta_x1)]);
    disp(['delta_x2 = ', num2str(delta_x2)]);

    % Display the updated values of x1 and x2 for this iteration
    disp(['Updated values:']);
    disp(['x1 = ', num2str(x1)]);
    disp(['x2 = ', num2str(x2)]);

    % Calculate percentage relative error for x1 and x2
    relative_error_x1 = abs((delta_x1 / x1) * 100);
    relative_error_x2 = abs((delta_x2 / x2) * 100);

    % Display the Jacobian inverse and percentage relative errors for x1 and x2
    disp(['Jacobian Inverse =']);
    disp(Jacobian_inverse);
    disp(['Percentage Relative Error for x1: ', num2str(relative_error_x1), '%']);
    disp(['Percentage Relative Error for x2: ', num2str(relative_error_x2), '%']);

    % Check for convergence
    if norm([delta_x1; delta_x2]) < tolerance
        break;
    end

    % Increment iteration counter
    iteration = iteration + 1;
end

% Print the final updated variables
disp(['Convergence achieved.']);
disp(['Final values:']);
disp(['x1 = ', num2str(x1)]);
disp(['x2 = ', num2str(x2)]);


% Plotting x1
figure;
plot(1:iteration, x1_values, 'b-o', 'LineWidth', 2);
grid on;
xlabel('Iteration');
ylabel('x1 Value');
title('Convergence of x1');

% Plotting x2
figure;
plot(1:iteration, x2_values, 'r-o', 'LineWidth', 2);
grid on;
xlabel('Iteration');
ylabel('x2 Value');
title('Convergence of x2');