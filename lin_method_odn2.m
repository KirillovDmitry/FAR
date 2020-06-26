% линейный метод, предложенный ќдносевцевым. ќцека амплитуды
% падабщей волны находитс€ как среднее значение от распределни€
% амплитуды пол€ по расрыву антенны. Ќет услови€ обрыва итераций --
% необходима априорна€ информаци€ о количестве волн.
function R = lin_method_odn2(Ampl, Phaza, X, m, lamda, tol, regl, ogr)

    % Ampl  - вектор амплитуды прин€того колебани€
    % Phaza - вектор фазы прин€того колебани€ (неразвернута€)
    % X     - вектор кооринат расположени€ вибраторов
    % m     - предполагаемое количество падающих волн (в алгоритме нет услови€ обрыва итераций)
    % lamda - длина волны
    % tol   - параметр в алгоритме развертки фазы unwrap


    eps = zeros(1, m);
    phi = zeros(1, m);
    ampl = zeros(1, m);
    P = zeros(1, m);

    for i=1:m

        Phaza = unwrap(Phaza, tol);
        [a, b] = lin_aproximate (Phaza, X);
        ampl(i) = find_amp(Ampl);

        phi(i) = b;
        eps(i) = asin(a*lamda/(2*pi));
        S = Ampl.*exp(1i.*Phaza) - ampl(i)*exp(1i*aproximate(2*pi*sin(eps(i))/lamda, phi(i), X));
        Ampl = abs( S );
        Phaza = angle( S );

    end
    N = m;
    IT = 0;
    [ampl_, eps_, phi_] = assort(ampl, eps, phi, m, regl);

    R = zeros(ogr,6);
    for i = 1:min(ogr,m)
        R(i,1) = eps_(i);
        R(i,2) = phi_(i);
        R(i,3) = ampl_(i);
        R(i,5) = P(i);
    end
    R(1,4) = N;
    R(1,6) = IT;