function  [eps_, phi_, ampl_, N, P, IT] = mp_one_method(Ampl, Phaza, X, m, lamda,eps_true, limit, N0,  H, h, regl)

    S = max(N0,m);
    P = zeros(1,m);
    phi    = -Inf(1,S);
    eps   = -Inf(1,S);
    ampl = -Inf(1,S);

    Y = Ampl(:).*exp(1i*Phaza(:));
    IT = 0;
    if m>0
        if norm(limit) == 0
            for i = 1:m
                [eps(i), ~,~, output] = fminsearch(@(xx)mp_one(xx, 1, X, lamda, Y), eps_true(i),  optimset(H, h(1), 'MaxIter', h(2)));
            end
            for i = 1:m
                ampl(i) = find_amp2(Y,eps(i), X, lamda);
            end
            IT =  output.iterations;
        else
            range1 = linspace(limit(2,1), limit(1, 2), N0+1);
            for i = 1:N0
                [eps(i),~, output] = my_fminsearch(@mp_one, 1, 0, X, lamda, Y, [range1(i) range1(i+1)], N0, H, h);
                ampl(i) = find_amp2(Y,eps(i), X, lamda);
                IT_ =  output.iterations;
                IT = IT + IT_;
            end
        end
    else
        m = 2;
        if norm(limit) == 0
            for i = 1:m
                [eps(i), ~,~, output] = fminsearch(@(xx)mp_one(xx, 1, X, lamda, Y), eps_true(i),  optimset(H, h(1), 'MaxIter', h(2)));
            end
            for i = 1:m
                ampl(i) = find_amp2(Y,eps(i), X, lamda);
            end
            IT =  output.iterations;
        else
            range1 = linspace(limit(2,1), limit(1, 2), N0+1);
            for i = 1:N0
                [eps(i),~, output] = my_fminsearch(@mp_one, 1, 0, X, lamda, Y, [range1(i) range1(i+1)], N0, H, h);
                ampl(i) = find_amp2(Y,eps(i), X, lamda);
                IT_ =  output.iterations;
                IT = IT + IT_;
            end
        end
        disp('mp_one_method: ветвь не отработана')
    end

    ampl_ = zeros(1, S);
    eps_ = zeros(1, S);
    phi_ = zeros(1, S);
    switch regl
        case 'ampl'
            [ampl,I] = sort(ampl);
            for i = 1:S
                eps_(i) = eps(I(length(ampl) - i +1));
                ampl_(i) = ampl(length(ampl) - i +1);
                phi_(i) = phi(I(length(ampl) - i +1));
            end
        case 'ygol'
            [eps,I] = sort(eps);
            for i = 1:S
                eps_(i) = eps(length(ampl) - i +1);
                ampl_(i) = ampl(I(length(ampl) - i +1));
                phi_(i) = phi(I(length(ampl) - i +1));
            end
    end

    eps_ = eps_(1:m);
    ampl_ = ampl_(1:m);
    phi_ = phi_(1:m);

    N =m;
    P(1)  =  sum( abs(Y).^2/2) / length(Ampl);
    
    
    if 0 % отладочная ветка
        NN = 100;
        II = linspace(-pi/5,pi/5,NN);
        G = zeros(1,NN);
        Y = Ampl(:).*exp(1i*Phaza(:));
        limit = [0 0.3; -0.3 0];
        for i=1:NN
            G(i) = mp_one(II(i), 1, X, lamda, Y);
        end
        
        figure(25); plot(II*180/pi,G);
        for mm = 1:m
            hold on;
            plot(eps_(mm)*180/pi, 1, '*');
        end
    end