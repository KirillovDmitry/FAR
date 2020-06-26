% метод максимального правдоподобия
function  [eps, phi, ampl] = mp_method(Ampl, Phaza, X, m, lamda, eps_true)

    phi  = zeros(1,2*ceil(m/2));
    eps  = zeros(1,2*ceil(m/2));
    ampl = zeros(1,2*ceil(m/2));

    Y1(:,1) = Ampl.*exp(1i*Phaza);
    Y2 = Y1';

    eps_true = eps_true +(rand(1,length(eps_true))-0.5).*eps_true*0.1;
    for i=1:ceil(m/2)
         
        [ eps(2*i-1: 2*i) ] = fminsearch(@(xx)mp(xx, [1,1], X, lamda, Y1, Y2),  eps_true(2*i-1: 2*i),  optimset('TolFun',1e-7, 'MaxIter',1000));
         if (nargin == 3)||(m > 2)
            [ampl(2*i-1: 2*i)] = fminsearch(@(xx)g(Y1, xx,eps(2*i-1: 2*i), X, lamda), [1, 0.8],  optimset('TolFun',1e-7, 'MaxIter',1000));
         else
             [ampl(2*i-1: 2*i)] = [1,1];
         end
   
         if (m > 2)
            S = Ampl.*exp(1i.*Phaza) - ampl(2*i-1)*exp(1i*2*pi*sin(eps(2*i-1))/lamda*X + phi(2*i-1));
            Ampl = abs( S );
            Phaza = angle( S );

            S = Ampl.*exp(1i.*Phaza) - ampl(2*i)*exp(1i*2*pi*sin(eps(2*i))/lamda*X + phi(2*i));
            Ampl = abs( S );
            Phaza = angle( S );
            Y1(:,1) = Ampl.*exp(1i*Phaza);
            Y2 = Y1';   
         end
         
    end
    
    
    function t = g(Z,a,eps, X, lamda)
        S1 =  exp(1i*2*pi*sin(eps(1))/lamda*X);
        S2 =  exp(1i*2*pi*sin(eps(2))/lamda*X);
        t = norm(Z(:) - a(1)*S1(:) - a(2)*S2(:));
        