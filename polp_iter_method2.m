% оптимизированный итерационный метод прямого-обратного линейного предсказания
function  R = polp_iter_method2(Ampl, Phaza, X, m, lamda, eps_true, limit,  N0, threshold, H, h, regl, ogr)

    Ampl_ = Ampl;
    Phaza_ = Phaza;

    if m>0
        eps = zeros(1,m);
        P = zeros(1,m);
        ampl = zeros(1,m);
        phi = zeros(1,m);
        eps_true = [eps_true 0 0 0 0 0];
        Y = Ampl(:).*exp(1i*Phaza(:));
        if norm(limit) == 0
            [eps(1), ~,~, output] = fminsearch(@(xx)polp(xx, 1, X, lamda, Y), eps_true(1),  optimset(H, h(1), 'MaxIter', h(2)));
        else
            [eps(1),~, output] = my_fminsearch(@polp, 1, phi, X, lamda, Y, [limit(1,1) limit(1,2)], N0, H, h);
        end
        IT =  output.iterations;
        ampl(1) = find_amp2(Y,eps(1), X, lamda);
        P(1) = sum( abs(Y).^2/2) / length(Ampl);

        for i=2:m
            S = Ampl.*exp(1i.*Phaza) - ampl(i-1)*exp(1i*2*pi*sin(eps(i-1))/lamda*X + phi(i-1));
            Ampl = abs( S );
            Phaza = angle( S );
            Y = Ampl(:).*exp(1i*Phaza(:));

            if norm(limit) == 0
                [eps(i), ~,~, output] = fminsearch(@(xx)polp(xx, 1, X, lamda, Y), eps_true(i),  optimset(H, h(1), 'MaxIter', h(2)));
            else
                [eps(i),~, output] = my_fminsearch(@polp, 1,  phi, X, lamda, Y, [limit(i,1) limit(i,2)], N0, H, h);
            end
            IT_ =  output.iterations;
            IT = IT + IT_;
            ampl(i) = find_amp2(Y,eps(i), X, lamda);
            P(i) = sum( abs(Y).^2/2) / length(Ampl);
        end
        N = m;
    else
        eps = zeros(1,1);
        ampl = zeros(1,1);
        phi = zeros(1,1);
        P = zeros(1,1);

        P1 = 10;  P2 = 20;     i = 1;    eps1 = 10;
        Y = Ampl(:).*exp(1i.*Phaza(:));
        P(1) = sum( abs(Y).^2/2) / length(Ampl);
        eps_true = [eps_true 0 0 0 0 0 0];
        if norm(limit) == 0
            [eps(1), ~,~, output] = fminsearch(@(xx)polp(xx, 1, X, lamda, Y), eps_true(i),  optimset(H, h(1), 'MaxIter', h(2)));
        else
            [eps(1),~, output] = my_fminsearch(@polp, 1, phi, X, lamda, Y, [limit(1,1) limit(1,2)], N0, H, h);
        end
        IT =  output.iterations;
        ampl(1) = find_amp2(Y,eps(1), X, lamda);

        while (P1>P(1)*threshold)  && (min(abs(eps(:)-eps1))>1*pi/180) && (abs(P1-P2)>P(1)*threshold) && (i<4)
            i = i +1;
            if i == 2
                P2 = P(i-1);
            else
                P2 = P1;
                phi(i-1) = 0;
                eps(i-1) = eps1;
                ampl(i-1) =    find_amp2(Y,eps1, X, lamda);
                P(i-1) = P1;
            end


            S = Ampl.*exp(1i.*Phaza) - ampl(i-1)*exp(1i*aproximate(2*pi*sin(eps(i-1))/lamda, phi(i-1), X));
            Ampl = abs( S );
            Phaza = angle( S );

            Y = Ampl(:).*exp(1i*Phaza(:));
            if norm(limit) == 0
                [eps1, ~,~, output] = fminsearch(@(xx)polp(xx,1, X, lamda, Y), eps_true(i),  optimset(H, h(1), 'MaxIter', h(2)));
            else
                [eps1,~, output] = my_fminsearch(@polp, 1, phi, X, lamda, Y, [limit(1,1) limit(1,2)], N0, H, h);
            end
            IT_ =  output.iterations;
            IT = IT + IT_;
            P1 = sum( abs(Y).^2/2) / length(Ampl);

        end
        N = i-1;

    end
    m = N;
    [ampl_, eps_, phi_] = assort(ampl, eps, phi,eps_true, m, regl);
    R = zeros(ogr,6);
    for i = 1:min(ogr,m)
        R(i,1) = eps_(i);
        R(i,2) = phi_(i);
        R(i,3) = ampl_(i);
        R(i,5) = P(i);
    end
    R(1,4) = N;
    R(1,6) = IT;

    if 0
        Y = Ampl_(:).*exp(1i*Phaza_(:));
        NN = 100;
        II = linspace(-pi/2,pi/2,NN);
        G = zeros(1,NN);
        amp0 = 1;
        for i=1:NN
            G(i) = polp(II(i), amp0, X, lamda, Y);
        end
        G = abs(G)./max(abs(G));
        figure(25);
        plot(II*180/pi,G);
        gg = max(G);
        hold on;
        plot(eps_(1)*180/pi, gg + 0.1*gg*(0.2+rand(1,1)), '.');
        hold on;
        plot(eps_(2)*180/pi, gg + 0.1*gg*(0.2+rand(1,1)), '.');
        hold on;
        plot(eps_true(1)*180/pi, gg + 0.2*gg, 'o');
        hold on;
        plot(eps_true(2)*180/pi, gg + 0.2*gg, 'o');
        %         figure(26);
        %         polar(exp(1i*II),G);
        %         hold on;
        %         polar(exp(1i*eps_(1)), 8, '*');
        %         hold on;
        %         polar(exp(1i*eps_(2)), 8, 'o');
    end

