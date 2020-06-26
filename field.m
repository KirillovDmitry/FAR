% рассчет распределения поля по раскрыву антенны с шумом и без него

function [Ampl_clear, Phaza_clear, Ampl_noise, Phaza_noise, Phaza_lin, P, sigma2] = field(eps, alpha, phi0, lamda, M, SNR, X)
	s = size(X);
	N = s(2);
	[Ampl, Phaza] = distrib(lamda, X(1), X(2)-X(1), N, M, eps, alpha, phi0);
	Ampl_clear = Ampl;                                         %   получили распределение поля по раскрыву антенны
	Phaza_clear = Phaza;                                       %   без шума
	Er = Phaza;
	Phaza_lin = unwrap(Phaza_clear);
	P = alpha(1)^2/2;
	if ischar(SNR)
		sigma2 = 0;
	else
		sigma2 =P/(10^(SNR/10));
	end
	noise1  =  sqrt(sigma2/2)*randn(1,N);                      %   добавляем шум в оба канала
	noise2  = sqrt(sigma2/2)*randn(1,N);

	noise_complex = (noise1+1i.*noise2);  
	pow_noise = sum(abs(noise_complex).^2)/N;  
	k = sigma2/pow_noise;
	noise1 = noise1*sqrt(k);
	noise2 = noise2*sqrt(k);

	a = real(Ampl_clear.*exp(1i.*Phaza_clear)) + noise1;
	b = imag(Ampl_clear.*exp(1i.*Phaza_clear)) + noise2;
	Ampl_noise = abs(a+1i.*b);                                %   получили распределение поля по раскрыву антенны
	Phaza_noise = angle(a+1i.*b);                             %   с шумом

	 
	for i = 1:N
		m = floor(abs(Phaza_noise(i) - Phaza_clear(i))/pi);
		Er(i) =  (Phaza_noise(i) - Phaza_clear(i)) - sign(Phaza_noise(i) - Phaza_clear(i))*m*2*pi;
	end

	Phaza_lin = Phaza_lin + Er;
