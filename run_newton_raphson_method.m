% Example usage
func = @(x) x^3 -10*x^2 +31.25*x -31.25;
derivative = @(x) 3*x^2-20*x+125/4;
initial_guess = 0;
tolerance = 0.0005;
true_value = 0.57940867;  % The true value you want to compare with

root = newton_raphson_method(tolerance, initial_guess, func, derivative, true_value);
if ~isnan(root)
    % Round off the root to 8 decimal digits
    rounded_root = round(root, 7);
    
    % Calculate the true error
    true_error = abs(true_value - rounded_root);
    
    % Calculate the percentage relative error
    percentage_relative_error = abs((rounded_root - initial_guess) / rounded_root) * 100;
    
    
end

function root = newton_raphson_method(tolerance, x0, func, derivative, true_value)
    % Implements the Newton-Raphson method to find a root of a given function.
    
    % Print initial value
    disp('Initial value:');
    disp(['x0 = ', num2str(x0, '%0.7f')]);
    disp(['f(x0) = ', num2str(func(x0), '%0.7f')]);  % Display f(x0)
    disp(['f_prime(x0) = ', num2str(derivative(x0), '%0.7f')]);  % Display f_prime(x0)
    
    iterations = 0;  % Start from 0 to match iteration count with display
    
    while true
        iterations = iterations + 1;
        
        % Compute the function value and derivative at x0
        f_x0 = func(x0);
        f_prime_x0 = derivative(x0);
        
        % Check if the derivative is close to zero
        if abs(f_prime_x0) < tolerance
            disp('The derivative is close to zero. Cannot continue the iteration.')
            root = NaN;
            return;
        end
        
        % Compute the next iteration value using Newton-Raphson formula
        x = x0 - f_x0 / f_prime_x0;
        
        % Round off the value of x to 8 decimal places
        rounded_x = round(x, 7);
        
        % Compute the function value and derivative at rounded_x
        f_x = func(rounded_x);
        f_prime_x = derivative(rounded_x);
        
        % Print current iteration and value
        disp(['Iteration ', num2str(iterations)]);
        disp(['x = ', num2str(rounded_x, '%0.7f')]);
        disp(['f(x) = ', num2str(f_x, '%0.7f')]);
        disp(['f_prime(x) = ', num2str(f_prime_x, '%0.7f')]);
        
        % Check if the value is already a root
        if f_x == 0
            root = rounded_x;
            return;
        end
        
        % Calculate the true error
        true_error = abs((true_value - rounded_x) / true_value) * 100;
        
        % Calculate the percentage relative error
        percentage_relative_error = abs((rounded_x - x0) / rounded_x) * 100;
        
        % Print true error and percentage relative error
        disp(['True error: ', num2str(round(true_error, 3))]);
        disp(['Percentage relative error: ', num2str(round(percentage_relative_error, 3))]);
        
        % Check if the percentage relative error is less than the tolerance
        if percentage_relative_error < tolerance
            break;
        end
        
        x0 = rounded_x;  % Update x0 for the next iteration
    end
    
    root = x0;
end
