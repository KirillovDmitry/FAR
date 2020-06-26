% линейный метод, предложенный ќдносевцевым. ќцека амплитуды
% падабщей волны находитс€ как среднее значение от распределни€
% амплитуды пол€ по расрыву антенны. Ќет услови€ обрыва итераций -
% необходима априорна€ информаци€ о количестве волн.
function [eps_, phi_, ampl_, N, P_, IT] = lin_method_odn(Ampl, Phaza, X, m, lamda, tol, regl)

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

    eps_ = zeros(1,N);
    P_ = P;
    ampl_ = zeros(1,N);
    phi_ = zeros(1,N);
    for i = 1:N
        ampl_(i) = ampl(i);
        eps_(i) = eps(i);
    end

    switch regl
        case 'ampl'
            [ampl,I] = sort(ampl);
            for i = 1:N
                if ampl(length(ampl) - i +1)~=0
                    ampl_(i) = ampl(length(ampl) - i +1);
                    eps_(i) = eps(I(length(ampl) - i +1));
                    phi_(i) = phi(I(length(ampl) - i +1));
                else
                    disp('ERROR: lin_method.m')
                end
            end
        case 'ygol'
            [eps,I] = sort(eps);
            for i = 1:N
                if ampl(length(ampl) - i +1)~=0
                    eps_(i) = eps(length(ampl) - i +1);
                    ampl_(i) = ampl(I(length(ampl) - i +1));
                    phi_(i) = phi(I(length(ampl) - i +1));
                else
                    disp('ERROR: lin_method.m')
                end
            end
    end
