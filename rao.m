function R = rao(a,b)
	f(a,b) = 0.5*a^2+2/b;
	R = diff(f, a);