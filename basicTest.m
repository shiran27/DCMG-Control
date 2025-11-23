clear; clc; close all;

%% -------------------------------------------------------------
% 1. TRUE SYSTEM (unknown to learner)
%% -------------------------------------------------------------
A = [0 1; -2 -3];
B = [0; 1];

n = size(A,1);
m = size(B,1);

% We will design state-feedback u = Kx
% true stabilizing K for checking:
K_true = [-2 -4];

%% -------------------------------------------------------------
% 2. SIMULATE CONTINUOUS TIME SYSTEM (ODE45)
%% -------------------------------------------------------------
Tfinal = 10;
dt_sim = 1e-3;
tspan = [0,1];

u_fun = @(t) 0.5*sin(2*t);   % persistently exciting input
w_fun = @(t) 0.01*randn(2,1); % process noise

x0 = [2;0];

f = @(t,x) A*x + B*u_fun(t) + w_fun(t);

[tsim, xsim] = ode45(f, tspan, x0);
tspan
xsim = xsim.';      % now 2 × N
usim = arrayfun(u_fun, tsim).';  % 1 × N

%% -------------------------------------------------------------
% 3. SAMPLE UNIFORMLY TO FORM DATA MATRICES
%% -------------------------------------------------------------
dt = 0.0001;  % data sampling interval
ts = tsim(1):dt:tsim(end);

x_s = interp1(tsim, xsim.', ts).';   % state samples 2×T
u_s = interp1(tsim, usim.', ts).';   % input samples 1×T

X  = x_s(:,1:end-1);       % x_k
Xp = x_s(:,2:end);         % x_{k+1}
U  = u_s(:,1:end-1);       % u_k

T = size(X,2);

%% -------------------------------------------------------------
% 4. DATA-DRIVEN CONTROLLER DESIGN
%    Using one-step relation:
%       Xp = A X + B U
%    Solve least-squares estimate:
%       [A B] = Xp * pinv([X; U])
%% -------------------------------------------------------------
Phi = [X; U];   % (n+m) × T
Theta = Xp * pinv(Phi);    % 2 × (2+1)

A
A_hat = Theta(:,1:n)
B
B_hat = Theta(:,n+1:end)

% Compute stabilizing LQR controller K = -R^{-1}B^T P
Q = eye(n);
R = 1;
[P_lqr,~,~] = care(A_hat, B_hat, Q, R);
K_dd = -R\(B_hat.'*P_lqr);

disp('True K:');
disp(K_true);
disp('Data-driven K:');
disp(K_dd);

%% -------------------------------------------------------------
% 5. VALIDATE CLOSED-LOOP BEHAVIOR
%% -------------------------------------------------------------
f_cl = @(t,x) A*x + B*K_dd*x + w_fun(t);

[ts2, xs2] = ode45(f_cl, [0 5], [2; 0]);

figure; hold on;
plot(ts2, xs2(:,1), 'LineWidth',2);
plot(ts2, xs2(:,2), 'LineWidth',2);
xlabel('t'); ylabel('states');
title('Closed-loop simulation with data-driven K');
legend('x_1','x_2');
grid on;
