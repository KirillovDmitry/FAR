% метод максимального правдоподобия одномерный (фактически измеритель мощности; сканирование лучом ДН)
function R = mp_one(eps, alpha, X, lamda,Y)
	S1 = alpha*exp(1i*2*pi*sin(eps)*X(:)/lamda);
	R = -sum(abs(S1'*Y)).^2/length(X);