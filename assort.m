    % функция производит сортировку входных массивов по убыванию угла,
    % амплитуды или близости углов к исходным в зависимости от параметра
    % regl.
	% m - длина массивов.
    function [Amp, Eps, Phi] = assort(amp, eps, phi, eps0, m, regl)

    amp_ = zeros(1, length(amp));
    phi_   = zeros(1, length(amp));
    eps_ = zeros(1, length(amp));

    Amp = zeros(1,m);
    Phi   = zeros(1, m);
    Eps  = zeros(1, m);

    switch regl
        case 'ampl'
             [amp, I] = unique(amp);
            for i = 1:length(amp)
                amp_(i) = amp(length(amp) - i +1);
                eps_(i)  = eps(I(length(amp) - i +1));
                phi_(i)   = phi(I(length(amp) - i +1));
            end
        case 'ygol'
            [eps,I] = unique(eps);
            for i = 1:length(eps)
                eps_(i)  = eps(length(eps) - i +1);
                amp_(i) = amp(I(length(eps) - i +1));
                phi_(i)   = phi(I(length(eps) - i +1));
            end
        case 'min'
            for i = 1:length(eps0)
                err = abs(eps - eps0(i));
                [~,B] = min(err);
                eps_(i)  = eps(B);
                amp_(i) = amp(B);
                phi_(i)   = phi(B);
                eps(B) = 1e10;
            end
        case 'norm'
            
          
    end

    m = min(length(eps), length(amp));
    j = 0;
    for i = 1:m
        if amp_(i)~=0
            j = j + 1;
            Eps(j) = eps_(i);
            Amp(j) = amp_(i);
            Phi(j) = phi_(i);
        end
    end

