% метод Кейпона
function R = kepon_one(eps, alpha, X, lamda,Y)
	S = alpha*exp(1i*2*pi*sin(eps)*X(:)/lamda);
	R = -1/(norm(S'*inv(Y)*S));
