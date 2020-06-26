% метод прямого-обратного линейного предсказания
function R = polp(eps, alpha, X, lamda,Y)
    L = 23;
    N = length(X);
    A = zeros(2*(N-L),L);
    h = zeros(2*(N-L),1);
    d2   = 1;
    for i = 1:2*(N-L)
        for j = 1:L
            if(i <= (N-L))
                A(i,j) = Y(i+L-j) ;
            else
                A(i,j) = conj(Y(i - N + L + j));
            end
        end

        if(i <= (N-L))
            h(i) = Y(L + i);
        else
            h(i) = conj(Y(i - (N - L)));
        end
    end
    AA = (A'*A)'*A';
    G = -AA*h;
    R = -abs(hh(sin(eps)*(2*pi*d2)/lamda, G))^2;
    
    function r = hh(phi,g)
        l = length(g);
        r = 0;
        for i = 1:l
            r = r+g(i)*(exp(1i*phi))^(-i);
        end
        r = r + 1;
