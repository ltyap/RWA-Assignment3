function coor = camberline(x, M, P)
%     M is the maximum camber divided by 100
%     P is the position of the maximum camber divided by 10
%     Zawsze u¿ywaj x sk³adaj¹cego siê ze 100 elementów do poprawnego
%     dzia³ania funkcji.
    x_front = x(1:10*P);
    x_back = x(10*P+1:end);
    M = M/100;
    P = P/10;
    yc_front = M/P^2*(2*P*x_front - x_front.^2);
    yc_back = M/(1-P)^2*(1-2*P+2*P*x_back-x_back.^2);
    coor = [yc_front, yc_back];
end