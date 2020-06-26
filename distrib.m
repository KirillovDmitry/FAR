% функция вычисления распределения поля (его аплитуды и фазы)
% по расскрыву фазированной антенной решетки 

function [Ampl, Phaza] = distrib(lamda, from, step, N, M, eps, alpha, phi0)
%	from - начальная координата Х = d/lamda
%	step - шаг по координате Х = d/lamda
%	N - количество точек рассчета по координате Х = d/lamda
%	М - количесво волн
%	eps(i) - угловая координата падения i-ой волны на ФАР
%	alpha(i) - амплитуда i-ой волны (нормируется на амплитуду первой волны)
%	Paza_turn - развернутая фаза на элементе ФАР
	
	I = 1;
    s = 0.0i;
    Ampl  = zeros(1,N);
    Phaza = zeros(1,N);
    
    for i=from:step:(from+step*(N-1))
         for j=1:M
            s = s + alpha(j)*exp(1.0000i*(phi0(j)+2*pi*i/lamda*sin(eps(j))));	
         end 
        Ampl(I) = abs(s);
        Phaza(I) =  angle(s);
        s = 0.0i; 
        I = I + 1;
    end

