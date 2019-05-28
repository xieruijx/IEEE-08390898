clear

load xi
[~, SampleNum] = size(xi_1);
[~, Index] = sort(rand(SampleNum, 1));
xi_1 = xi_1(:, Index);

save('xirandom', 'xi_1', 'xi_2_Set');