% функция вычисления ошибки оценки угловой координаты
    function [err var] = error_(eps, eps_method, m, alpha, regl)

    err = zeros(1,m);
    var = zeros(1,m);
    eps_ = zeros(1,m);
    j = 0;
    for i = 1:m
        if alpha(i)~=0
            j = j + 1;
            eps_(j) = eps(i);
        end
    end
    
    switch regl
        case 'ampl'
            [alpha,I] = sort(alpha);
            I = flipdim(I,2);
            for i = 1:length(eps_method)
                eps_(i) = eps(I(i));
            end
            eps = eps_;
        case 'ygol'
            eps = eps_;
            eps = sort(eps);
            eps = flipdim(eps,2);
    end
    if m==1
        err(1) =  ((eps_method(1)) - eps(1) );
        var(1) = ((eps_method(1)) - eps(1) )^2;
    else
        
        for  i=1:m
            err(i) = eps_method(i)  - eps(i);
            var(i) = (eps_method(i)  - eps(i) )^2;
        end
    end

