        % ������� ��������� �������� ����������� �� ������ ��, 
        % ������� �����-�� ������, �� ��� �������������� �������� � 
        % �� ����� (��������� ��������� ��������� � ����������� �� ������
        % ���������, �� ������ ������� ��������� � ���� �������
        % ����� ������ ���������)
		
function R = mp(eps, alpha, X, lamda,Y1, Y2)
      
      	% eps   - ������-������� ����� ������� ����� (���������)
        % alpha - ������-������� �������� ����� (�� ��������� �� ���������)
        % X     - ������ �������� ������������ ����������
        % lamda - ����� �����
        % Y1    - ������-������� ��������� ���������
        % Y2    - ������-������ ��������� ���������, ����������� �� ������ ������������ Y1
        
    S1 = alpha(1).*exp(1i*2*pi*sin(eps(1)).*X/lamda);
    S2 = alpha(2).*exp(1i*2*pi*sin(eps(2)).*X/lamda);			
    F1(:,1) = S1;
    F1(:,2) = S2;
    F2 = F1';
    R = -norm(Y2*(F1/(F2*F1)*F2)*Y1);