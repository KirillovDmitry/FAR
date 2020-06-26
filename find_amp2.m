% алгоритм поиска амплитуды падающего сигнала 
% из одномерного метода максимального правдоподобия
function ampl = find_amp2(A, eps, X, lamda)
	S = exp(1i*2*pi*sin(eps)*X(:)/lamda);
	ampl = norm(S'*A/length(X));