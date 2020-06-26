function R=psi_function( eps, X, phi0, alpha, lamda)

	phi1 = phi0(2) + 2*pi*X/lamda*sin(eps(1));
	phi2 = phi0(1) + 2*pi*X/lamda*sin(eps(2));
	d_phi = phi2 - phi1;
	alp = alpha(2)/alpha(1);
	a = 1 + alp*cos(d_phi);
	b = alp*sin(d_phi);

	c = angle( 1j*(a.*sin(phi1)+b.*cos(phi1))+ (a.*cos(phi1)-b.*sin(phi1)));
	R = c;