% одномерный итерационный метод максимального правдоподобия
% в данном методе производится оценка количества падающих лучей по 
% результатам измерения остаточной мощности сигнала
function  [eps_, phi_, ampl_, N, P, IT] = mp_one_iter_method(Ampl, Phaza, X, m, lamda, eps_true, limit,  N0, threshold, H, h, regl)
   
    if m>0
        eps = zeros(1,m);
        P = zeros(1,m);
        ampl = zeros(1,m);
        phi = zeros(1,m);
        Ampl_ = Ampl;
        Phaza_ = Phaza;
        eps_true = [eps_true 0 0 0 0 0];
        Y = Ampl(:).*exp(1i*Phaza(:));
        if norm(limit) == 0           
            [eps(1), ~,~, output] = fminsearch(@(xx)mp_one(xx, 1, X, lamda, Y), eps_true(1),  optimset(H, h(1), 'MaxIter', h(2)));
        else
            [eps(1),~, output] = my_fminsearch(@mp_one, 1, phi, X, lamda, Y, [limit(1,1) limit(1,2)], N0, H, h);
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
                [eps(i), ~,~, output] = fminsearch(@(xx)mp_one(xx, 1, X, lamda, Y), eps_true(i),  optimset(H, h(1), 'MaxIter', h(2)));
            else
                [eps(i),~, output] = my_fminsearch(@mp_one, 1,  phi, X, lamda, Y, [limit(i,1) limit(i,2)], N0, H, h);
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
            [eps(1), ~,~, output] = fminsearch(@(xx)mp_one(xx, 1, X, lamda, Y), eps_true(i),  optimset(H, h(1), 'MaxIter', h(2)));
        else
            [eps(1),~, output] = my_fminsearch(@mp_one, 1, phi, X, lamda, Y, [limit(1,1) limit(1,2)], N0, H, h);
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
                ampl(i-1) = find_amp2(Y,eps1, X, lamda);
                P(i-1) = P1;
            end

            S = Ampl.*exp(1i.*Phaza) - ampl(i-1)*exp(1i*aproximate(2*pi*sin(eps(i-1))/lamda, phi(i-1), X));
            Ampl = abs( S );
            Phaza = angle( S );

            Y = Ampl(:).*exp(1i*Phaza(:));
            if norm(limit) == 0
                [eps1, ~,~, output] = fminsearch(@(xx)mp_one(xx, 1, X, lamda, Y), eps_true(i),  optimset(H, h(1), 'MaxIter', h(2)));
            else
                [eps1,~, output] = my_fminsearch(@mp_one, 1, phi, X, lamda, Y, [limit(1,1) limit(1,2)], N0, H, h);
            end

            IT_ =  output.iterations;
            IT = IT + IT_;
            P1 = sum( abs(Y).^2/2) / length(Ampl);

        end
        N = i-1;
        m = N;
    end
    
    eps_ = zeros(1,m);
    ampl_ = zeros(1,m);
    phi_ = zeros(1,m);
    
    switch regl
        case 'ampl'
            [ampl,I] = sort(ampl);
            for i = 1:m
                if ampl(length(ampl) - i +1)~=0
                    eps_(i) = eps(I(length(ampl) - i +1));
                    ampl_(i) = ampl(length(ampl) - i +1);
                    phi_(i) = phi(I(length(ampl) - i +1));
                else
                    disp('ERROR: mp_one_iter_method.m')
                end
            end
        case 'ygol'
            [eps,I] = sort(eps);
            for i = 1:m
                eps_(i) = eps(length(ampl) - i +1);
                ampl_(i) = ampl(I(length(ampl) - i +1));
                phi_(i) = phi(I(length(ampl) - i +1));
            end
    end
    
    eps_ = eps_(1:m);
    ampl_ = ampl_(1:m);
    phi_ = phi_(1:m);
    
    if 0
        NN = 100;
        II = linspace(-pi/5,pi/5,NN);
        G = zeros(1,NN);
        Y = Ampl_(:).*exp(1i*Phaza_(:));
        limit = [0 0.3; -0.3 0];
        for i=1:NN
            G(i) = mp_one(II(i), 1, X, lamda, Y);
        end
        
        a2 = max(abs(G));
        a1 = min(abs(G));
        G = (abs(G)-a1)/(a2-a1);
        figure(25); plot(II*180/pi,G);
        for mm = 1:m
            hold on;
            plot(eps_(mm)*180/pi, 1.1, '*');
        end
    end