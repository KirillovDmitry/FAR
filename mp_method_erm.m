% метод МП, описанный в методичке Ермолаева
function  [eps_, phi_, ampl_, N, P_, IT] = mp_method_erm(Ampl, Phaza, X, lamda, eps_true, limit,N0, H, h, regl)
    

    m = 2; % оригиналный метод двумерный
    phi    = zeros(1,2*ceil(m/2));
    eps   = ones(1,2*ceil(m/2));
    ampl = zeros(1,2*ceil(m/2));
    P = zeros(1,2*ceil(m/2));
    Z = Ampl(:).*exp(1i*Phaza(:));
    eps_true = [eps_true 0 0 0 0 0 0];
    Ampl_ = Ampl;
    Phaza_ = Phaza;
        
    if norm(limit) == 0
        [eps(1: 2), ~,~, output] = fminsearch(@(xx)mp_erm(xx, [1 1], X, lamda, Z),  eps_true(1: 2),  optimset(H, h(1), 'MaxIter', h(2)));
    else
        [eps(1: 2),~, output] = my_fminsearch(@mp_erm, [1 1], phi, X, lamda, Z, limit, N0, H, h);
    end
    IT =  output.iterations;
    if (nargout >= 3)
        [ampl(1: 2)] = fminsearch(@(xx)g(Z, xx,eps(1: 2), X, lamda), [1, 1],  optimset(H, h(1), 'MaxIter', h(2)));
        P(1) = sum( abs(Z).^2/2) / length(Ampl);
        S = Ampl.*exp(1i.*Phaza) - ampl(1)*exp(1i*2*pi*sin(eps(1))/lamda*X + phi(1));
        Ampl = abs( S );
        Phaza = angle( S );
        Z = Ampl(:).*exp(1i*Phaza(:));
        P(2) = sum( abs(Z).^2/2) / length(Ampl);
    end
    
    N = 2;
    eps_ = zeros(1,length(ampl));
    P_ = P;
    ampl_ = zeros(1,length(ampl));
    phi_ = zeros(1,length(ampl));
    switch regl
        case 'ampl'
            [ampl,I] = sort(ampl);
            for i = 1:length(ampl)
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
            for i = 1:length(ampl)
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
        II = linspace(-pi/12,pi/12,NN);
        G = zeros(NN,NN);
        Y = Ampl_(:).*exp(1i*Phaza_(:));
        for i=1:NN
            for j=1:NN
                G(i,j) = mp_erm([II(i), II(j)],[1 1], X, lamda, Y);
            end
        end
        figure(25); 
        contour(II*180/pi,II*180/pi,G, 100);
        hold on;
        plot3(eps_(1)*180/pi, eps_(2)*180/pi,8, '*');
        hold on;
        plot3(eps_true(1)*180/pi, eps_true(2)*180/pi,8, 'o');
        
        
    end
       
    function t = g(Z,a,eps, X, lamda)
    S1 =  exp(1i*2*pi*sin(eps(1))/lamda*X);
    S2 =  exp(1i*2*pi*sin(eps(2))/lamda*X);
    t = norm(Z(:) - a(1)*S1(:) - a(2)*S2(:));
    
    
    