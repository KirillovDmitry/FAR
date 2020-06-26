function R=psi45(eps, X, lamda, Z) % phi0 не используется

	alpha = eps(3);
	d_phi = eps(4);
	eps  = eps(1:2);
	s = exp(1i*(2*pi*X/lamda*sin(eps(1)))) + alpha(1)*exp(1i*(d_phi(1)+2*pi*X/lamda*sin(eps(2))));	  
	A = abs(s);
	P =  angle(s);
   
	c = Z - A(:).*exp(1i.*P(:));
	R = norm(c);
