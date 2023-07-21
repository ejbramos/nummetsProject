function [V, theta] = newtonRaphson3Bus(Ybus)

    % Step 1: Define the power system data
    P = [0; 1.5; -1.0]; % Active power injections (P) in per unit (p.u.)
    Q = [0; 0.5; -0.5]; % Reactive power injections (Q) in per unit (p.u.)
    V = [1.0; 1.0; 1.0]; % Initial voltage magnitude guesses (p.u.)
    theta = [0; 0; 0]; % Initial voltage phase angle guesses (radians)

    % Step 2: Set tolerance and maximum iterations
    tol = 1e-6;
    maxIter = 50;

    % Step 3: Implement the Newton-Raphson method
    iter = 0;
    while iter < maxIter
        iter = iter + 1;

        % Step 3a: Calculate power injections and admittances for each bus
        [Pcalc, Qcalc, G, B] = calculatePowerInjectionsAndAdmittances(V, theta, Ybus);

        % Step 3b: Calculate mismatch vector
        Pmismatch = P - Pcalc;
        Qmismatch = Q - Qcalc;
        mismatch = [Pmismatch; Qmismatch];

        % Step 3c: Calculate Jacobian matrix
        J = calculateJacobian(G, B, V, theta, Ybus);

        % Step 3d: Update voltage magnitudes and phase angles
        update = J \ mismatch;
        dTheta = update(1:3);
        dV = update(4:6);
        theta = theta + dTheta;
        V = V + dV;

        % Step 3e: Check for convergence
        if max(abs(mismatch)) < tol
            break;
        end
    end

    % Step 4: Display the results
    fprintf('Converged in %d iterations.\n', iter);
    fprintf('Bus 1: V = %.4f p.u., theta = %.4f radians.\n', V(1), theta(1));
    fprintf('Bus 2: V = %.4f p.u., theta = %.4f radians.\n', V(2), theta(2));
    fprintf('Bus 3: V = %.4f p.u., theta = %.4f radians.\n', V(3), theta(3));
end

function [Pcalc, Qcalc, G, B] = calculatePowerInjectionsAndAdmittances(V, theta, Ybus)
    % Step 3a: Calculate power injections and admittances for each bus
    G = real(Ybus);
    B = imag(Ybus);

    S = V .* (cos(theta) + 1i * sin(theta)); % Complex bus voltages
    I = Ybus * S; % Complex bus currents

    Pcalc = real(S .* conj(I)); % Real power injections
    Qcalc = imag(S .* conj(I)); % Reactive power injections
end

function J = calculateJacobian(G, B, V, theta, Ybus)
    % Step 3c: Calculate Jacobian matrix
    n = length(V);
    J = zeros(2 * n);

    for i = 1:n
        for j = 1:n
            if i == j
                J(i, j) = -Qcalc(i) - V(i)^2 * B(i, i);
                J(i, j + n) = Pcalc(i) - V(i)^2 * G(i, i);
                J(i + n, j) = Pcalc(i) + V(i) * V(j) * (G(i, j) * sin(theta(i) - theta(j)) - B(i, j) * cos(theta(i) - theta(j)));
                J(i + n, j + n) = -Pcalc(i) + V(i) * V(j) * (G(i, j) * cos(theta(i) - theta(j)) + B(i, j) * sin(theta(i) - theta(j)));
            else
                J(i, j) = V(i) * V(j) * (G(i, j) * sin(theta(i) - theta(j)) - B(i, j) * cos(theta(i) - theta(j)));
                J(i, j + n) = V(i) * V(j) * (G(i, j) * cos(theta(i) - theta(j)) + B(i, j) * sin(theta(i) - theta(j)));
                J(i + n, j) = V(i) * V(j) * (G(i, j) * cos(theta(i) - theta(j)) + B(i, j) * sin(theta(i) - theta(j)));
                J(i + n, j + n) = -V(i) * V(j) * (G(i, j) * sin(theta(i) - theta(j)) - B(i, j) * cos(theta(i) - theta(j)));
            end
        end
    end
end

% Given Ybus matrix
Ybus = [53.8516*exp(-1i*68.1986), 22.3606*exp(1i*116.5651), 31.6228*exp(1i*108.4349);
        22.3606*exp(1i*116.5651), 58.1378*exp(-1i*63.4349), 35.7771*exp(1i*116.5651);
        31.6228*exp(1i*108.4349), 35.7771*exp(1i*116.5651), 67.2309*exp(-1i*67.2490)];

% Call the load flow analysis function with the given Ybus
[V_final, theta_final] = newtonRaphson3Bus(Ybus);
