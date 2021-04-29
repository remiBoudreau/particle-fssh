% Fewest Switches Surface Hopping Simple Avoided Crossing Model by Tully,
% 1990, JCP
%
% Jean-Michel Boudreau June 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Values presented in atomic units.
 
clc
clear
global A B C D c1 c2 x n m PES c1_dot c2_dot V11 i
% Constants for Potentials - Change as Desired
A = 0.01; 
B = 1.6;
C = 0.005;
D = 1.0;

% Parameters - Change as Desired
n_sim = 100; % Set number of trajectories
x_0 = -5; % Initial position
dt = 1; % Time step 
m = 2000; % Mass
p_vector = [6];
disp(length(p_vector)) % For determining progress (end of program); See line 
% Check for probability of transmission dependance on average momentum
for p_in = 1:length(p_vector) % Average momentum
    p_0 = p_vector(p_in);
    % Set up initial values for all 'i' trajectories
    x = zeros(n_sim, 10000);
    p = zeros(n_sim, 10000);
    c1 = zeros(n_sim, 10000);
    c2 = zeros(n_sim, 10000);
    tm1 = 0;
    tm2 = 0;
    ref = 0;
    sigma_x = 20/p_0; % Standard Deviation of x 
    sigma_p = p_0/40; % Standard Deviation of p considering uncertainty principle sigma_x*sigma_p >= h_bar/2
    % Normalized distribution of x and p values for all 'i' trajectories
    x(:,1) = randn(n_sim, 1);
    x(:,1) = x(:,1)/std(x(:, 1))*sigma_x; 
    x(:,1) = x(:,1) - mean(x(:, 1)) + x_0;
    p(:,1) = randn(n_sim, 1);
    p(:,1) = p(:,1)/std(p(:, 1))*sigma_p; 
    p(:,1) = p(:,1) - mean(p(:, 1)) + p_0;
    c1(:, 1) = 1;
    % Set up trajectories
    for i = 1:n_sim
        PES = 1; % PES under consideration; PES = 1 is G.S., PES = 2 is E.S.
        % Set-up for ith trajectory 
        c1_dot = 0;
        c2_dot = 0;
        y_n = [x(i, 1); p(i, 1); c1(i, 1); c2(i, 1)];
        n = 1;
        if x(i, 1) <= 0
            V11 = -exp(-D*x(i, n)^2)*(C^2 + A^2*exp(2*D*x(i, n)^2) - 2*A^2*exp(B*x(i, n))*exp(2*D*x(i, n)^2) + A^2*exp(2*B*x(i, n))*exp(2*D*x(i, n)^2))^(1/2); % Adiabatic
        elseif x(i, 1) > 0
            V11 = -exp(-B*x(i, n))*exp(-D*x(i, n)^2)*(C^2*exp(2*B*x(i, n)) + A^2*exp(2*D*x(i, n)^2) - 2*A^2*exp(B*x(i, n))*exp(2*D*x(i, n)^2) + A^2*exp(2*B*x(i, n))*exp(2*D*x(i, n)^2))^(1/2); % Adiabatic
        end
        
        if p(i, 1)^2 > 2*m*(-0.005 - V11) % Specific condition for paramters A B C D for more efficiency. If KE < (max(V11) - V11(x_0)) then no transmission can occur 
            if (x(i, 1) < 0 && p(i, 1) > 0) || (x(i, 1) > 0 && p(i, 1) < 0)
                % Surface Hopping Algorithm
                while (x(i, n) <  10) && (x(i, n) >  -10) 
                         % Surface Hop for nth iteration for use in n + 1 iteration 
                          if PES == 1
                              a = abs(c1(i, n))^2;
                              a_dot = conj(c2_dot)*c2(i, n) + conj(c2(i, n))*c2_dot;
                              parm = dt*a_dot/a;
                              if parm > rand && p(i, n)^2 > abs(4*m*V11);  % Surface hopping criterion for 1 -> 2
                                  PES = 2; % Surface hop
                                  y_n(2) = sqrt(p(i, n)^2 + 4*m*V11);
                              end
                        elseif PES == 2
                            a = abs(c2(i, n))^2;
                            a_dot = conj(c1_dot)*c1(i, n) + conj(c1(i, n))*c1_dot;
                            parm = dt*a_dot/a;
                            if parm > rand % Surface hopping criterion for 2 -> 1
                                PES = 1; % Surface hop
                                y_n(2) = sqrt(p(i, n)^2 - 4*m*V11);
                            end
                          end

                        % Runge-Kutta Integration
                        k_1 = dt * RK_SH(y_n);
                        k_2 = dt * RK_SH(y_n + k_1/2);
                        k_3 = dt * RK_SH(y_n + k_2/2);
                        k_4 = dt * RK_SH(y_n + k_3);
                        x(i, n + 1) = y_n(1) + k_1(1)/6 + k_2(1)/3 + k_3(1)/3 + k_4(1)/6;
                        p(i, n + 1) = y_n(2) + k_1(2)/6 + k_2(2)/3 + k_3(2)/3 + k_4(2)/6;
                        c1(i, n + 1) = y_n(3) + k_1(3)/6 + k_2(3)/3 + k_3(3)/3 + k_4(3)/6;
                        c2(i, n + 1) = y_n(4) + k_1(4)/6 + k_2(4)/3 + k_3(4)/3 + k_4(4)/6;
                        y_n = [x(i, n + 1); p(i, n + 1); c1(i, n + 1); c2(i, n + 1)];
                        n = n + 1;
    %                     % Conservation of Energy Check: Potential Energy for Total Energy below 
    %                     if PES == 1
    %                         PE(i, n) = V11;
    %                     elseif PES == 2
    %                         PE(i, n) = -V11;
    %                     end
                end
            end
        end
        % Count number of trajectories that transmit
        if (x(i, n) > 0 && x(i, 1) < 0) || (x(i, n) < 0 && x(i, 1) > 0)
            if PES == 1
                tm1 = tm1 + 1;
            elseif PES == 2
                tm2 = tm2 + 1;
            end
        else
            ref = ref + 1;
        end
    end
    
    % Calculate probabilities of transmission and reflection
    prob_tm1(1, p_in) = tm1/n_sim;
    prob_tm2(1, p_in) = tm2/n_sim;
    prob_ref(1, p_in) = ref/n_sim;
    disp(p_in); % For determining progress (by present index of average momentums)
%     Conservation of energy and state population check
%     t = 1:10000;
%     c1_2 = abs(c1).^2;
%     c2_2 = abs(c2).^2;
%     pop = c1_2 + c2_2;
%     KE = p.^2/2;
%     KE(:, end) = [];
%     plot(t, pop)
%     for i = 1:n_sim
%         figure(i)
%         plot(t, KE(i, :), t, PE(i, :))
%         figure(i + n_sim)
%         plot(t, c1_2(i, :), t, c2_2(i, :), t, pop(i, :))
%     end
end

% Graph %%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
plot(p_vector, prob_tm1, 'b')
hold on
plot(p_vector, prob_tm2, 'r')
hold off
figure(2)
plot(p_vector, prob_ref, 'm')
