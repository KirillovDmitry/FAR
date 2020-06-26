% �������� �����, ������������ ������������. ����� ���������
% �������� ����� ��������� ��� ������� �������� �� ������������
% ��������� ���� �� ������� �������. ��� ������� ������ �������� --
% ���������� ��������� ���������� � ���������� ����.
function R = lin_method_odn2(Ampl, Phaza, X, m, lamda, tol, regl, ogr)

    % Ampl  - ������ ��������� ��������� ���������
    % Phaza - ������ ���� ��������� ��������� (�������������)
    % X     - ������ �������� ������������ ����������
    % m     - �������������� ���������� �������� ���� (� ��������� ��� ������� ������ ��������)
    % lamda - ����� �����
    % tol   - �������� � ��������� ��������� ���� unwrap


    eps = zeros(1, m);
    phi = zeros(1, m);
    ampl = zeros(1, m);
    P = zeros(1, m);

    for i=1:m

        Phaza = unwrap(Phaza, tol);
        [a, b] = lin_aproximate (Phaza, X);
        ampl(i) = find_amp(Ampl);

        phi(i) = b;
        eps(i) = asin(a*lamda/(2*pi));
        S = Ampl.*exp(1i.*Phaza) - ampl(i)*exp(1i*aproximate(2*pi*sin(eps(i))/lamda, phi(i), X));
        Ampl = abs( S );
        Phaza = angle( S );

    end
    N = m;
    IT = 0;
    [ampl_, eps_, phi_] = assort(ampl, eps, phi, m, regl);

    R = zeros(ogr,6);
    for i = 1:min(ogr,m)
        R(i,1) = eps_(i);
        R(i,2) = phi_(i);
        R(i,3) = ampl_(i);
        R(i,5) = P(i);
    end
    R(1,4) = N;
    R(1,6) = IT;