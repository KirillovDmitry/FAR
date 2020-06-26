% линейный метод  Односевцева, усовершенствованный. Оцека амплитуды
% падабщей волны находится из одномерного алгоритма МП.
% Условие обрыва итерации - падание остаточной мощности ниже порога.

function [eps_, phi_, ampl_, N, P_, IT] = lin_method(Ampl, Phaza, X, m, lamda, tol, threshold, regl, lin)
    % Ampl      - вектор амплитуды принятого колебания
    % Phaza     - вектор фазы принятого колебания (неразвернутая)
    % X         - вектор кооринат расположения вибраторов
    % m         - предполагаемое количество падающих волн (в алгоритме нет условия обрыва итераций)
    % lamda     - длина волны
    % tol       - параметр в алгоритме развертки фазы unwrap
    % threshold - порог по мощности для окончания итерационного процесса
    IT = 0;
    if m>0 % если знаем сколько должно прийти волн
        eps = zeros(1,m);
        ampl = zeros(1,m);
        phi = zeros(1,m);
        P =  zeros(1,m);

        for i=1:m

            if lin == 1
                Phaza = unwrap(Phaza, tol);
            else 
                if i ~= 1
                    Phaza = unwrap(Phaza, tol);
                end
            end
            Z = Ampl(:).*exp(1i.*Phaza(:));
            [a, b] = lin_aproximate(Phaza, X);
            phi(i) = b;
            eps(i) = asin(a*lamda/(2*pi));
            ampl(i) = find_amp2(Z,eps(i), X, lamda);
            P(i) = sum( abs(Z).^2/2) / length(Ampl);
            S = Ampl.*exp(1i.*Phaza) - ampl(i)*exp(1i*aproximate(2*pi*sin(eps(i))/lamda, phi(i), X));

            Ampl = abs( S );
            Phaza = angle( S );
        end
        N = m;
    else % если не знаем сколько должно прийти волн, итерационный алгоритм

        eps = zeros(1,1);
        ampl = zeros(1,1);
        phi = zeros(1,1);
        P =  zeros(1,1);
        Z = Ampl(:).*exp(1i.*Phaza(:));
        P1 = 10;  P2 = 20;     i = 1;    eps1 = 10;
        P(1) = sum( abs(Z).^2/2) / length(Ampl);
        if lin == 1
            Phaza = unwrap(Phaza, tol);            
        end
        [a, b] = lin_aproximate(Phaza, X);
        phi(1) = b;
        eps(1) = asin(a*lamda/(2*pi));
        ampl(1) = find_amp2(Z,eps(1), X, lamda);

        while (P1>P(1)*threshold)  && (min(abs(eps(:)-eps1))>1*pi/180) && (abs(P1-P2)>P(1)*threshold)
            i = i +1;

            if i == 2
                P2 = P;
            else
                P2 = P1;
                phi(i-1) = b;
                eps(i-1) = eps1;
                ampl(i-1) = find_amp2(Z,eps1, X, lamda);
                P(i-1) = P2;
            end

            S = Ampl.*exp(1i.*Phaza) - ampl(i-1)*exp(1i*aproximate(2*pi*sin(eps(i-1))/lamda, phi(i-1), X));
            Ampl = abs( S );
            Phaza = angle( S );

            Phaza = unwrap(Phaza, tol);
            Z = Ampl(:).*exp(1i.*Phaza(:));
            [a, b] = lin_aproximate(Phaza, X);
            eps1 = asin(a*lamda/(2*pi));
            P1 = sum( abs(Z).^2/2) / length(Ampl);
        end
        N = i-1;
    end


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
    
 