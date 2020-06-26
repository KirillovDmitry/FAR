% метод Кейпона
function R = kepon_one_method(KM, Y, X, m, lamda, eps_true, limit, N0, threshold, H, h, regl, ogr)

    S = max(N0,m);
    P = zeros(1, S);
    phi    = zeros(1, S);
    eps   = zeros(1, S);
    ampl = zeros(1, S);
    P(1)  =  sum( abs(Y).^2/2) / length(KM);

    IT = 0;
    if m>0
        if norm(limit) == 0
            for i = 1:m
                [eps(i), ~,~, output]  = fminsearch(@(xx)kepon_one(xx, 1,  X, lamda, KM), eps_true(i),  optimset(H, h(1), 'MaxIter', h(2)));
                IT_ =  output.iterations;
                IT = IT + IT_;
            end
            for i = 1:m
                ampl(i) = find_amp2(Y,eps(i), X, lamda);
            end
        else
            range1 = linspace(limit(2,1), limit(1, 2), N0+1);
            for i = 1:N0
                [eps(i),~, output] = my_fminsearch(@kepon_one, 1, 0, X, lamda, KM, [range1(i) range1(i+1)], N0, H, h);
                ampl(i) = find_amp2(Y,eps(i), X, lamda);
                IT_ =  output.iterations;
                IT = IT + IT_;
            end
        end
    else
        range1 = linspace(limit(2,1), limit(1, 2), N0+1);
        for i = 1:N0
            [eps(i),~, output] = my_fminsearch(@kepon_one, 1, 0, X, lamda, KM, [range1(i) range1(i+1)], N0, H, h);
            ampl(i) = find_amp2(Y,eps(i), X, lamda);
            if ampl(i) < h(3)*sqrt(2*P(1))
                ampl(i) = 0;
                eps(i) = 0;
            end
            IT_ =  output.iterations;
            IT = IT + IT_;
        end
    end

    [ampl_, eps_, phi_] = assort(ampl, eps, phi, m, regl);
    N = length(ampl_);
    
    R = zeros(ogr,6);
    for i = 1:min(ogr,length(ampl_))
        R(i,1) = eps_(i);
        R(i,2) = phi_(i);
        R(i,3) = ampl_(i);
        R(i,5) = P(i);
    end
    R(1,4) = N;
    R(1,6) = IT;
    
    
    if 0
        NN = 100;
        II = linspace(-pi/5,pi/5,NN);
        G = zeros(1,NN);
        limit = [0 0.3; -0.3 0];
        for i=1:NN
            G(i) = kepon_one(II(i), 1, X, lamda, KM);
        end

        g = max(abs(G));
        figure(25); plot(II*180/pi,G/g);
        for mm = 1:m
            hold on;
            plot(eps_(mm)*180/pi, 1, '*');
        end
    end
