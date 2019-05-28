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
alpha_h = calculateAlpha(alpha, d_0)

%% 优化
% 变量
x_1 = sdpvar(1, 1);
x_2 = sdpvar(1, 1);

% w_0 = sdpvar(1, 1);
p_1 = sdpvar(T, q);
p_2 = sdpvar(T, q);
p_3 = sdpvar(T, q);

f = sdpvar(q, 1);
g = sdpvar(q, 1);
gamma = sdpvar(1, 1);

% 目标函数
Objective = Pr_1 * x_1 + Pr_2 * x_2  ;

% 约束
Constraint1 = [x_1 >= 0, x_2 >= 0,  w_0 >= 0];

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
    Constraint7 = [Constraint7, xi_2(:, k) - p_1(:, k) - p_3(:, k) - f(k) <= 0];
end
Constraint8 = [g >= f - gamma * ones(q, 1)];

Constraint9 = [g >= 0];

Constraint10 = [gamma + 1 / q / alpha_h * sum(g) <= 0];

Constraints = [Constraint1, Constraint2, Constraint3, Constraint4, Constraint5, Constraint6, Constraint7, Constraint8, Constraint9, Constraint10];

% 计算
sol = optimize(Constraints, Objective);
% Analyze error flags
if sol.problem == 0
 % Extract and display value
 solution_x_1 = value(x_1);
 solution_x_2 = value(x_2);

 [solution_x_1, solution_x_2]
%  SOC_0 = value(w_0) / solution_x_2
 solution_p_1 = value(p_1);
 solution_p_2 = value(p_2);
 solution_p_3 = value(p_3);

 solution_Objective = value(Objective)
 solution_gamma = value(gamma)
 solution_f = value(f);
 solution_g = value(g);
 maxf = max(solution_f);
 minf = min(solution_f);
 maxg = max(solution_g);
 ming = min(solution_g);
 [maxf, minf, maxg, ming]
 solution_gamma + 1 / q / alpha_h * sum(solution_g)
 rate = 1 / q * sum(double(solution_f > 0))
%  solution_p_1 + solution_p_2 - xi_1 * solution_x_1
%  xi_2 - solution_p_1 - solution_p_3
else
 display('Hmm, something went wrong!');
 sol.info
 yalmiperror(sol.problem)
end
toc
x1 = solution_x_1;
x2 = solution_x_2;