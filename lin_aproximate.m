% метод линейной аппроксимации (метод наименьших квадратов)
function [a, b, err] = lin_aproximate (Y, X)

    M_x = 0; M_y = 0; M_xy = 0; M_xx = 0;
    err = 0;
    s = size(X);

    for  i = 1: 1: s(2)
        M_x = M_x + X(i);
        M_y = M_y + Y(i);
        M_xx = M_xx + X(i)^2;
        M_xy = M_xy + X(i)*Y(i);
    end

    M_x    =   M_x/ s(2);
    M_y    =   M_y/ s(2);
    M_xx  =   M_xx/ s(2);
    M_xy  =   M_xy/ s(2);

    a = (M_xy-M_x*M_y)/(M_xx-M_x*M_x);
    b = M_y - a*M_x;
  
    for i = 1:  s(2)
        err = err +  (a*X(i) + b - Y(i))^2;
    end
    err = err/(s(2));
    