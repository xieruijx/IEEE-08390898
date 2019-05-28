clear

load xirandom
q = 200;
alpha = 0.05;
d_0 = 0.02;

[~, ~, StationNum] = size(xi_2_Set);

x1_Set = zeros(StationNum, 1);
x2_Set = zeros(StationNum, 1);

xi_1 = xi_1(:, 1: q);
for i = 1: StationNum
    xi_2 = xi_2_Set(:, :, i);
    [x1, x2, alpha_h, rate] = getPlanningResult(q, alpha, d_0, xi_1, xi_2)
    x1_Set(i) = x1;
    x2_Set(i) = x2;
end

Capacity = [x1_Set, x2_Set]

save('Capacity', 'Capacity');