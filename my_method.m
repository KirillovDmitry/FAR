% несколько оптимизированный линейный метод односевцева
function [eps_, phi_, ampl_, N, P_, IT] = my_method(Ampl, Phaza, X, m, lamda, tol, threshold, regl)
	IT = 0;
	N = 1;
	Phaza = unwrap(Phaza, tol);
	Z = Ampl(:).*exp(1i.*Phaza(:));
	[a, b] = lin_aproximate(Phaza, X);
	phi_(1) = b;
	eps_(1) = asin(a*lamda/(2*pi));
	ampl_(1) = find_amp2(Z,eps, X, lamda);
	P_ = sum( abs(Z).^2/2) / length(Ampl);
