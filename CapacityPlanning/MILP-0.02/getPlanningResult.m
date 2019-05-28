function [x1, x2, alpha_h, rate] = getPlanningResult(q, alpha, d_0, xi_1, xi_2)
tic
T = 24;
Pr_1 = 1;
Pr_2 = 1;
w_l = 0.1;
w_h = 0.9;
p_2m = 0.2;
p_3m = 0.2;
eta_1 = 0.9;
eta_2 = 0.9;
w_0 = 0.3;
delta_t = 1;
M = 2500;
alpha_h = calculateAlpha(alpha, d_0)

% 变量
x_1 = sdpvar(1, 1);
x_2 = sdpvar(1, 1);

% w_0 = sdpvar(1, 1);
p_1 = sdpvar(T, q);
p_2 = sdpvar(T, q);
p_3 = sdpvar(T, q);

z = binvar(q, 1);

% 目标函数
Objective = Pr_1 * x_1 + Pr_2 * x_2 ;

% 约束
Constraint1 = [x_1 >= 0, x_2 >= 0, w_0 >= 0];

Constraint2 = [];
for i = 1: q
    Constraint2 = [Constraint2, p_1(:, i) >= 0, p_2m * x_2 >= p_2(:, i) >= 0, p_3m * x_2 >= p_3(:, i) >= 0];
end

allonestril = tril(ones(T));
Constraint3 = [];
for i = 1: q
    Constraint3 = [Constraint3, w_l * x_2 * ones(T, 1) <= w_0 * x_2 + allonestril * (eta_1 * p_2(:, i) - 1 / eta_2 * p_3(:, i)) * delta_t <= w_h * x_2 * ones(T, 1)];
end

Constraint4 = [];
for i = 1: q
    Constraint4 = [Constraint4, sum(eta_1 * p_2(:, i) - 1 / eta_2 * p_3(:, i)) == 0 ];
end

Constraint5 = [];

Constraint6 = [];
for k = 1: q
    Constraint6 = [Constraint6, p_1(:, k) + p_2(:, k)  - xi_1(:, k) * x_1 <= 0];
end

Constraint7 = [];
for k = 1: q
    Constraint7 = [Constraint7, xi_2(:, k) - p_1(:, k) - p_3(:, k) <= M * z(k)];
end

Constraint8 = [sum(z) <= q * alpha_h];

Constraints = [Constraint1, Constraint2, Constraint3, Constraint4, Constraint5, Constraint6, Constraint7, Constraint8];

% 计算
% options = sdpsettings('verbose',1,'solver','mosek');
options = sdpsettings('verbose',1);
sol = optimize(Constraints, Objective, options);
% Analyze error flags
if sol.problem == 0
 % Extract and display value
 solution_x_1 = value(x_1);
 solution_x_2 = value(x_2);

 [solution_x_1, solution_x_2]
 solution_w_0 = value(w_0)
 solution_p_1 = value(p_1);
 solution_p_2 = value(p_2);
 solution_p_3 = value(p_3);

 solution_Objective = value(Objective)
 solution_z = value(z);
 rate = 1 / q * sum(solution_z)
%  size(xi_2)
%  size(solution_p_1)
%  size(solution_p_3)
 shouldbenegative = max(max(xi_2(:, 1: q) - solution_p_1 - solution_p_3 - M))
%  solution_p_1 + solution_p_2 - xi_1 * solution_x_1
%  xi_2 - solution_p_1 - solution_p_3
else
 display('Hmm, something went wrong!');
 sol.info
 yalmiperror(sol.problem)
end
toc
x1 = value(x_1);
x2 = value(x_2);