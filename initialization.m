	%%  инициализация основных массивов для хранения результатов моделирования
    for temp = 1:ogr
        for met = 1:12
            evalin('base', ['Amp_method',int2str(met),'_', int2str(temp),' = zeros(Mean,N_z);']);
            evalin('base', ['Angle_method',int2str(met),'_', int2str(temp),' = zeros(Mean,N_z);']);
            evalin('base', ['Power_method',int2str(met),'_', int2str(temp),' = zeros(Mean,N_z);']);
            evalin('base', ['Err_method',int2str(met),'_', int2str(temp),' = zeros(Mean,N_z);']);
            evalin('base', ['Sko_method',int2str(met),'_', int2str(temp),' = zeros(Mean,N_z);']);
            evalin('base', ['SKO_method',int2str(met),'_', int2str(temp),' = zeros(1,N_z);']);
        end
    end
	% основной массив для хранения результатов моделирования
    Matrix2 = zeros(6, Met, ogr, Mean, N_z); % измеряемая величина, номер метода, номер луча, усреднение, количество итераций по варьируемому параметру)
    null = 'NULL';

    for met = 1:12
        evalin('base', ['Iteration',int2str(met), '= 0;']);
        evalin('base', ['Iter',int2str(met), '= zeros(Mean,N_z);']);
        evalin('base', ['IT',int2str(met), '= 0;']);
        evalin('base', ['NumberSignal',int2str(met), '= 0;']);
        evalin('base', ['NSignal',int2str(met), '=  zeros(Mean,N_z);']);
        evalin('base', ['NS',int2str(met), '= 0;']);
        evalin('base', ['P',int2str(met), '= 0;']);
        evalin('base', ['Err_method',int2str(met), '= zeros(1,ogr);']);
        evalin('base', ['method',int2str(met), ' = null;']);
    end
	% основной массив для хранения результатов моделирования
    Matrix1 = zeros(9, Met,  Mean, N_z); % измеряемая величина, номер метода, усреднение, количество итераций по варьируемому параметру)

    for met = 1:12
        evalin('base', 'eps_method = zeros(met,m);');
        evalin('base', 'phi_method = zeros(met,m);');
        evalin('base', 'ampl_method = zeros(met,m);');
    end

    MethodBool = zeros(1,12);
    MethodStr  = ['method01'; 'method02';  'method03'; 'method04'; 'method05';  'method06'; 'method07';  'method08'; 'method09';  'method10'; 'method11';  'method12'];
    NOISE = zeros(Mean,N_z);
    Q = zeros(Mean,N_z);
    KM = zeros(N,N);
    NumberSignal = 0;
    T1= zeros(Mean,N_z);
    T2= zeros(Mean,N_z);
    T3= zeros(Mean,N_z);
