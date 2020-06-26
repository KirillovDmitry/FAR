% фазовый метод оценки параметров падающей волны (четырехмерный)
function [eps_, phi_, ampl_,N, P_, IT] = psi4_method(Ampl, Phaza, X, m, lamda, eps_true, limit, N0, threshold, H, h, regl)
    Ampl_ = Ampl;
    Phaza_ = Phaza;
    amp0 = [1 0.95];
    if m>0
        phi    = zeros(1,2*ceil(m/2));
        P = zeros(1,2*ceil(m/2));
        eps_true = [eps_true 0 0 0 0 0 0];
        Z = Ampl(:).*exp(1i*Phaza(:));
        P(1) = sum( abs(Z).^2/2) / length(Ampl);

        if norm(limit) == 0
            [temp, ~,~, output] = fminsearch(@(xx)psi4(xx, phi(1:2), X, lamda, Z), [eps_true(1:2) amp0],  optimset(H, h(1), 'MaxIter', h(2)));
            eps = temp(1:2);
            ampl = temp(3:4);
        else
            [temp,~, output] = my_fminsearch(@psi4,  amp0, phi, X, lamda, Z, limit, N0, H, h);
            eps = temp(1:2);
            ampl = temp(3:4);
        end
        IT =  output.iterations;
        
        S = Ampl.*exp(1i.*Phaza) - ampl(1)*exp(1i*2*pi*sin(eps(1))/lamda*X + phi(1));
        Ampl = abs( S );
        Phaza = angle( S );
        Z = Ampl(:).*exp(1i*Phaza(:));
        P(2) = sum( abs(Z).^2/2) / length(Ampl);
        S = Ampl.*exp(1i.*Phaza) - ampl(2)*exp(1i*2*pi*sin(eps(2))/lamda*X + phi(2));
        Ampl = abs( S );
        Phaza = angle( S );
        Z = Ampl(:).*exp(1i*Phaza(:));


        for i=2:ceil(m/2)
            P(2*i-1) = sum( abs(Z).^2/2) / length(Ampl);
            S = Ampl.*exp(1i.*Phaza) - ampl(2*i-1-2)*exp(1i*2*pi*sin(eps(2*i-1-2))/lamda*X + phi(2*i-1-2));
            Ampl = abs( S );
            Phaza = angle( S );
            Z = Ampl(:).*exp(1i*Phaza(:));
            P(2*i) = sum( abs(Z).^2/2) / length(Ampl);
            S = Ampl.*exp(1i.*Phaza) - ampl(2*i-2)*exp(1i*2*pi*sin(eps(2*i-2))/lamda*X + phi(2*i-2));
            Ampl = abs( S );
            Phaza = angle( S );
            Z = Ampl(:).*exp(1i*Phaza(:));
            
            if norm(limit) == 0
                [temp, ~,~, output] = fminsearch(@(xx)psi4(xx, phi(2*i-1: 2*i), X, lamda, Z), [eps_true(2*i-1: 2*i) amp0],  optimset(H, h(1), 'MaxIter', h(2)));
                eps(2*i-1: 2*i) = temp(1:2);
                ampl(2*i-1: 2*i) = temp(3:4);
            else
                [temp,~, output] = my_fminsearch(@psi4,  amp0, phi, X, lamda, Z, limit, N0, H, h);
                eps(2*i-1: 2*i) = temp(1:2);
                ampl(2*i-1: 2*i) = temp(3:4);
            end
            IT_ =  output.iterations;
            IT = IT + IT_;

        end
        N = m;

    else

        eps_true = [eps_true 0 0 0 0 0 0];
        P = zeros(1,2);

        P12 = 10;  P21 = 20; P11 = 30;  P22 = 40;    i = 1;    eps1 = 10;   eps2 = 10;
        Z = Ampl(:).*exp(1i.*Phaza(:)); ampl1 = 1;  ampl2 = 1;
        P(1) = sum( abs(Z).^2/2) / length(Ampl);
        phi = [0 0 0 0 0 0 0 0];

        if norm(limit) == 0
            [temp, ~,~, output] = fminsearch(@(xx)psi4(xx, phi(1:2), X, lamda, Z), [eps_true(1:2) amp0],  optimset(H, h(1), 'MaxIter', h(2)));
            eps = temp(1:2);
            ampl = temp(3:4);
        else
            [temp,~, output] = my_fminsearch(@psi4,  amp0, phi, X, lamda, Z, limit, N0, H, h);
            eps = temp(1:2);
            ampl = temp(3:4);
        end
        IT =  output.iterations;

        while (P11>P(1)*threshold) && (P12>P(1)*threshold) && (min(abs(eps(:)-eps1))>1*pi/180)  && (min(abs(eps(:)-eps2))>1*pi/180) && (abs(P21-P22)>P(1)*threshold) && (abs(P11-P12)>P(1)*threshold) && (i<4)
            i = i +1;
            if i == 2
                P21 = P(2*i-2);
                P22 = P(2*i-3);
            else
                P21 = P11;
                P22 = P12;
                [eps(2*i-1-2: 2*i-2) ] =  [eps1 eps2];
                [ampl(2*i-1-2: 2*i-2)] = [ampl1 ampl2];
                [phi(2*i-1-2: 2*i-2)] = [0 0];
                P(2*i-1-2) = P11;
                P(2*i-2) = P12;
            end

            S = Ampl.*exp(1i.*Phaza) - ampl(2*i-1-2)*exp(1i*2*pi*sin(eps(2*i-1-2))/lamda*X + phi(2*i-1-2));
            Ampl = abs( S );
            Phaza = angle( S );
            Z = Ampl(:).*exp(1i.*Phaza(:));
            P11 = sum( abs(Z).^2/2) / length(Ampl);

            S = Ampl.*exp(1i.*Phaza) - ampl(2*i-1-1)*exp(1i*2*pi*sin(eps(2*i-1-1))/lamda*X + phi(2*i-1-1));
            Ampl = abs( S );
            Phaza = angle( S );
            Z = Ampl(:).*exp(1i.*Phaza(:));

            if norm(limit) == 0
                [temp, ~,~, output] = fminsearch(@(xx)psi4(xx, phi(1:2), X, lamda, Z), [eps_true(1:2) amp0],  optimset(H, h(1), 'MaxIter', h(2)));
                e = temp(1:2);
                a = temp(3:4);
                ampl1 = a(1);
                ampl2 = a(2);
            else
                [temp,~, output] = my_fminsearch(@psi4,  amp0, phi, X, lamda, Z, limit, N0, H, h);
                e = temp(1:2);
                ampl1 = e(3);
                ampl2 = e(4);
            end
            IT_ =  output.iterations;
            IT = IT + IT_;

            eps2  = e(2);
            eps1 = e(1);
            P12 = sum( abs(Z).^2/2) / length(Ampl);

        end
        N = (i -1)*2;
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
                    disp('ERROR: mp_method_erm.m')
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
                    disp('ERROR: mp_method_erm.m')
                end
            end
    end

    if 0
        NN = 100;
        II = linspace(-pi/5,pi/5,NN);
        G = zeros(NN,NN);
        Y = Ampl_(:).*exp(1i*Phaza_(:));
        for i=1:NN
            for j=1:NN
                G(i,j) = psi4([II(i), II(j), 1, 0.8], [0 0], X, lamda, Y);
            end
        end
        [~, t] = min(abs(II-eps_(2)));
        m2 = max(abs(G(:, t)));
        m1 = min(abs(G(:, t)));
        a2 = max(max(abs(G)));
        a1 = min(min(abs(G)));
        T = (abs(G(:, t))-m1)/(m2-m1);
        G = (abs(G)-a1)/(a2-a1);
                figure(23); mesh(II*180/pi,II*180/pi, G);
        hold on;

        plot3(eps_(1)*180/pi, eps_(2)*180/pi, 1.1, 'ro');
        hold on;
        plot3(eps_true(1)*180/pi, eps_true(2)*180/pi, 1.1, '*');
        figure(24); hold on; plot(II*180/pi, T, 'r');
        hold on; plot(eps_(1)*180/pi,1.1, 'ro');
        
%         figure(25); mesh(II*180/pi,II*180/pi,G);
%         hold on;
%         plot3(eps_(1)*180/pi, eps_(2)*180/pi,8, '*');
    end
