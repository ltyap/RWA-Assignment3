function eta = parabCamber(x, eps)
eta = zeros(1, length(x));
    for i = 1:length(x)
        eta(i) = 4*eps*x(i)*(1-x(i));
    end
end