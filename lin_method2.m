% �������� �����  �����������, �������������������. ����� ���������
% �������� ����� ��������� �� ����������� ��������� ��.
function R = lin_method2(Ampl, Phaza, X, m, lamda, eps0, tol, threshold, regl, lin, ogr)
	% Ampl       --   ������ ��������� ��������� ���������
	% Phaza      --   ������ ���� ��������� ��������� (�������������)
	% X          --   ������ �������� ������������ ����������
	% m          --   �������������� ���������� �������� ���� (� ��������� ��� ������� ������ ��������)
	% lamda      --   ����� �����
	% tol        --   �������� � ��������� ��������� ���� unwrap
	% threshold  --   ����� �� �������� ��� ��������� ������������� ��������

    IT = 0;
    if m>0 % ���� ����� ������� ������ ������ ����
        eps = zeros(1,m);
        ampl = zeros(1,m);
        phi = zeros(1,m);
        P =  zeros(1,m);

        for i=1:m
            if norm(lin) == 1
                Phaza = my_unwrap(Phaza, tol);
            else
                if i == 1
                    Phaza = lin;
                else
                    Phaza = my_unwrap(Phaza, tol);
                end
            end
            Z = Ampl(:).*exp(1i.*Phaza(:));
            [a, b] = lin_aproximate(Phaza, X);
            phi(i) = b;
            eps(i) = asin(a*lamda/(2*pi));
            ampl(i) = find_amp2(Z,eps(i), X, lamda);
            P(i) = sum( abs(Z).^2/2) / length(Ampl);
            S = Ampl.*exp(1i.*Phaza) - ampl(i)*exp(1i*aproximate(2*pi*sin(eps(i))/lamda, phi(i), X));

            Ampl = abs( S );
            Phaza = angle( S );
        end
    else % ���� �� ����� ������� ������ ������ ����, ������������ ��������
    end
    N = m;
    [ampl_, eps_, phi_] = assort(ampl, eps, phi, eps0, m, regl);
    R = zeros(ogr,6);
    for i = 1:min(ogr,m)
        R(i,1) = eps_(i);
        R(i,2) = phi_(i);
        R(i,3) = ampl_(i);
        R(i,5) = P(i);
    end
    R(1,4) = N;
    R(1,6) = IT;



