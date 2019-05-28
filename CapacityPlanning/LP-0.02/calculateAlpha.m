function alpha_h = calculateAlpha(alpha, d_0)
    fun = @(x)(exp(-d_0) * power(x, 1-alpha) - 1) / (x - 1);
    [~, alpha_c] = fmincon(fun, alpha, [], [], [], [], 0, 1);
    alpha_h = 1 - alpha_c;