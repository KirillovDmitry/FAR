% ������� ���������� ������������� ���� (��� �������� � ����)
% �� ��������� ������������ �������� ������� 

function [Ampl, Phaza] = distrib(lamda, from, step, N, M, eps, alpha, phi0)
%	from - ��������� ���������� � = d/lamda
%	step - ��� �� ���������� � = d/lamda
%	N - ���������� ����� �������� �� ���������� � = d/lamda
%	� - ��������� ����
%	eps(i) - ������� ���������� ������� i-�� ����� �� ���
%	alpha(i) - ��������� i-�� ����� (����������� �� ��������� ������ �����)
%	Paza_turn - ����������� ���� �� �������� ���
	
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

