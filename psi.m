function R=psi(eps, phi0, alpha, X, lamda, Z)
    s = alpha(1)*exp(1i*(phi0(1)+2*pi*X/lamda*sin(eps(1)))) + alpha(2)*exp(1i*(phi0(2)+2*pi*X/lamda*sin(eps(2))));	  
    R = sum(abs(Z - s(:)).^2);
