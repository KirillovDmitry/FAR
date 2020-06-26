 function R=psi4(eps, phi0, X, lamda, Z)
	alpha = eps(3:4);
	eps  = eps(1:2);
	s = alpha(1)*exp(1i*(phi0(1)+2*pi*X/lamda*sin(eps(1)))) + alpha(2)*exp(1i*(phi0(2)+2*pi*X/lamda*sin(eps(2))));	  
	A = abs(s);
	P =  angle(s);
	c = Z - A(:).*exp(1i.*P(:));
	R = norm(c);
