	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% скрипт предназначен для моделирования алгоритмов приема переотраженного зондирующего
	% сигнала фазированной антенной решеткой. проводится сравнительный анализ по оценке 
	% эффективности нескольких методов оценки угловой координаты цели, подсвеченной антенной.
	% ввод основных параметров, а также вывод результатов моделирования, реализован с
	% помощью стандартых средств Matlab'a по разработки графического интерфейса.
	
	%% 
    clear all; close all; clc;
    %matlabpool open local 4 % запуск нескольких потоков
    %%  задаем характеристики распределения поля по раскрыву антенны

	% задание основных переменных моделирования
    M = 2;                                           % количество сигналов, падающих на ФАР
    m = 2;                                           % предполагаемое количество сигналов, падающих на ФАР
    alpha1 = 1.0;     phi1 = 0;     eps1 = 5;        % задаем фазы, углы падения (в градусах) и амплитуды
    alpha2 = 0.5;     phi2 = 0;     eps2 = -5;       % падающих сигналов.
    alpha3 = 0.0;     phi3 = 160;   eps3 = -3;       % амплитуды нормированны на амплитуду первой волны
    alpha4 = 0.0;     phi4 = -53;   eps4 = 0;

    threshold = 0.1;
    lamda = 1.61;                                    % отношение  d/lamda=0.66 (по ТЗ)
    d = 1;                                           % при варьировании величины х/lamda геометрическое
                                                     % местоположение вибраторов соответствует её
                                                     % целым значениям
													  
    SNR = 12;  %'none'		                         % задаем отношение сигнал шум по мощности;
													 % мощность первой волны принимаем равной 0.5
    tol = 1*pi;

    phi0   = [phi1,phi2, phi3, phi4 0 0 0 0]*pi/180; % вектора фаз, амплитуд и углов падения волн
    alpha  = [alpha1,alpha2, alpha3, alpha4 0 0 0 0];
    eps    = [eps1,eps2, eps3, eps4 0 0 0 0]*pi/180;


    from = 0;                                        % пределы моделирования по координате Х = x/lamda ФАР
    step = 1;
    N = 30;                                          % количество точек рассчета по координате Х = d/lamda
    X = from:step:(from+step*(N-1));

    from_z = 1;                                      % пределы моделирования по координате Z
    to_z = 30;
    N_z = 40;                                        % количество точек моделирования по координате Z
    Z = linspace(from_z, to_z, N_z);
    method = 'SNR' ; % MaxIter, SNR, a, tol EPS1, EPS2, EPS_sim, EpsReflected, PHI2, precision; % углы в градусах
    Mean = 500;                                      % количество усреднений по шуму
    Met = 12;                                        % количечтво методов
    ogr = 8;                                         % ограничение на колличество выводимых волн
    H =  'TolFun';                                   % TolFun, TolX
    h1 =  [1e-5 2000 0.3];                           % h(3) - порог по мощности
    N0 = 7;                                          % разбиение предела на N0 интервалов
   
	
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
	
	%% вспомогательное окно для выбора алгоритма оценивания координат цели с помощью ФАР
    MethodBool = GetMethod();
	
	% инициализация переменных моделирования
    zz = 0;
    h = waitbar(0,'Please wait...');            TStart = tic;
    dT = 0;
    T0 = datestr(clock);
    screenSize = get(0,'ScreenSize');
    pointsPerPixel = 1;
    width = 500 * pointsPerPixel;
    height = 100 * pointsPerPixel;
    pos = [screenSize(3)/2-width/2 screenSize(4)/2+height/1 width height];
    handles.WinBar = figure('MenuBar', 'None',...
        'BusyAction', 'queue',...
        'WindowStyle','normal',...
        'Interruptible', 'off', ...
        'DockControls', 'off', ...
        'Name', 'Time',...
        'NumberTitle', 'Off', ...
        'Visible','off',...
        'Position',pos);
    Txt = uicontrol('Style', 'Text',...
        'String','Hello, world!',   'Position', [50, 15, 400, 70]);
    str = ['начало моделирования: '  T0];
    set(Txt, 'String', str);
    set(handles.WinBar, 'Visible',    'on');
	% главный цикл вариации параметра
    for z = Z                                                     
        zz = zz+1;
        for mm = 1:Mean                  % цикл усреднения
            Ttemp = tic;
			% в зависимости от выбранного метода формируем массив параметра моделирования
            switch method 
                case 'SNR'
                    SNR = Z(zz);
                case 'tol'
                    tol = Z(zz)*pi/180;
                case 'a'
                    alpha(2) = Z(zz);
                case 'EPS1'
                    eps(1) = Z(zz)*pi/180;
                case 'EPS2'
                    eps(2) = Z(zz)*pi/180;
                case 'EPS_sim'
                    eps(1) = Z(zz)*pi/180;
                    eps(2) = -Z(zz)*pi/180;
                case 'EPS_reflected'
                    eps(1) = Z(zz)*pi/180;
                    eps(2) = -Z(zz)*pi/180+1.5*pi/180;
                case 'precision'
                    h1(1) = Z(zz);
                case 'MaxIter'
                    Z = round(Z);
                    h1(2) = Z(zz);
                otherwise
                    disp('ошибка выбора метода')
                    pause;
            end

            %% вычисляем распределение поля по раскрыву антенны
            [Ampl_clear, Phaza_clear, Ampl_noise, Phaza_noise, Phaza_lin, Pow, sigma2] = field(eps, alpha, phi0, lamda, M, SNR, X);
            Phaza_turn_clear    =   unwrap(Phaza_clear);
            Phaza_turn_noise   =   unwrap(Phaza_noise, tol);
            UnwrapError = norm( (Phaza_turn_noise-Phaza_turn_clear).^2 );

            %% вычисляем сгенерированный уровень шума в сигнале
            temp = sum( abs(Ampl_noise.*exp(1i*Phaza_noise) - Ampl_clear.*exp(1i*Phaza_clear)).^2) / length(Ampl_clear);
            NOISE(mm,zz) = 10*log10( (Pow) / temp );

            if strcmp('SNR', method)
                SnrLevel = abs(from_z - to_z)/N_z*0.4;
            else
                SnrLevel = SNR*0.05;
            end

            if ischar(SNR)
                [Ampl_clear, Phaza_clear, Ampl_noise, Phaza_noise, Pow, sigma2] = field(eps, alpha, phi0, lamda, M, SNR, X);
                temp = sum( abs(Ampl_noise.*exp(1i*Phaza_noise) - Ampl_clear.*exp(1i*Phaza_clear)).^2) / length(Ampl_clear);
                NOISE(mm,zz) = 10*log10( (Pow) / temp );
            else
                while( abs(NOISE(mm,zz) - SNR) > SnrLevel)
                    [Ampl_clear, Phaza_clear, Ampl_noise, Phaza_noise, Pow, sigma2] = field(eps, alpha, phi0, lamda, M, SNR, X);
                    temp = sum( abs(Ampl_noise.*exp(1i*Phaza_noise) - Ampl_clear.*exp(1i*Phaza_clear)).^2) / length(Ampl_clear);
                    NOISE(mm,zz) = 10*log10( (Pow) / temp );
                end
            end
            Q(mm,zz) = sqrt(Pow/temp);
            Phaza_turn_clear   =   unwrap(Phaza_clear);       % развернутая фаза входного сигнала без шума
            Phaza_turn_noise   =   unwrap(Phaza_noise, tol);  % развернутая фаза входного сигнала с шумом

            %% определяем координаты сигналов различными методами
            Ampl_temp       =  Ampl_noise;
            Phaza_temp     =  Phaza_noise;

            eps0 = eps;
			limit = [-0*pi/180 0*pi/180; -0*pi/180 0*pi/180];
            regl = 'ygol'; %'ygol' 'ampl'

            if MethodBool(1) == 1;
                regl_temp = 'ampl'; 
                N0_ = 15;
                if 1
                    [eps_method1, phi_method1, ampl_method1, NumberSignal1, P1, Iteration1] =...
                        mp_one_method(Ampl_temp, Phaza_temp, X, m, lamda, eps0, limit, N0, H, h1, regl_temp);
                    method1 = 'MpOneMethod';
                else
					% отладочные варианты
                    [eps_method1, phi_method1, ampl_method1, NumberSignal1, P1, Iteration1] =...
                        mp_one_method(Ampl_temp, Phaza_temp, X, m, lamda, eps, limit, N0, H, h1, regl_temp);
                    method1 = 'A';
                    [eps_method2, phi_method2, ampl_method2, NumberSignal2, P2, Iteration2] =...
                        mp_one_method(Ampl_temp, Phaza_temp, X, m, lamda, eps0, limit, N0, H, h1, regl_temp);
                    method2 = 'B';
                    limit = [0*pi/180 30*pi/180; -30*pi/180 0*pi/180; 0*pi/180 30*pi/180; -30*pi/180 0*pi/180];
                    [eps_method3, phi_method3, ampl_method3, NumberSignal3, P3, Iteration3] =...
                        mp_one_method(Ampl_temp, Phaza_temp, X, m, lamda, eps0, limit, N0, H, h1, regl_temp);
                    method3 = 'C';
                end
            end

            if MethodBool(2) == 1;
				[eps_method2, phi_method2, ampl_method2, NumberSignal2, P2, Iteration2] =...
					mp_one_iter_method(Ampl_temp, Phaza_temp, X, m, lamda,eps0, limit, N0, threshold, H, h1, regl);
				method2 = 'MpOneIterMethod';
            end

            if MethodBool(3) == 1;
				[eps_method3, phi_method3, ampl_method3, NumberSignal3, P3, Iteration3] =...
					mp_method_erm(Ampl_temp, Phaza_temp, X, lamda,eps0,limit, N0, H, h1, regl);
				method3 = 'MpMethodErm';
            end

            if MethodBool(4) == 1;
                if 1
                    [eps_method4, phi_method4, ampl_method4, NumberSignal4, P4, Iteration4] =...
                        mp_iter_method_erm(Ampl_temp, Phaza_temp, X, m, lamda, eps0, limit, N0, threshold, H, [1e-3 1000],regl);
                    method4 = 'MpIterMethodErm';
                else
                    [eps_method1, phi_method1, ampl_method1, NumberSignal1, P1, Iteration1] =...
                        mp_iter_method_erm(Ampl_temp, Phaza_temp, X, m, lamda, eps,limit, N0, threshold, H, h1,regl);
                    method1 = 'A';
                    eps_method1
                    [eps_method2, phi_method2, ampl_method2, NumberSignal2, P2, Iteration2] =...
                        mp_iter_method_erm(Ampl_temp, Phaza_temp, X, m, lamda, eps0,limit, N0, threshold, H, h1,regl);
                    method2 = 'B';
                    limit = [0*pi/180 30*pi/180; -30*pi/180 0*pi/180; 0*pi/180 30*pi/180; -30*pi/180 0*pi/180];
                    [eps_method3, phi_method3, ampl_method3, NumberSignal3, P3, Iteration3] =...
                        mp_iter_method_erm(Ampl_temp, Phaza_temp, X, m, lamda, eps0,limit, N0, threshold, H, h1,regl);
                    method3 = 'C';
                end
            end

            if MethodBool(5) == 1;
                [eps_method5, phi_method5, ampl_method5, NumberSignal5, P5, Iteration5] =...
                    lin_method_odn(Ampl_temp, Phaza_temp, X, m, lamda, tol, regl);
                method5 = 'LinMethodOdn';

            end

            if MethodBool(6) == 1;
                [eps_method6, phi_method6, ampl_method6, NumberSignal6, P6, Iteration6] =...
                    lin_method(Ampl_temp, Phaza_temp, Phaza_lin, X, m, lamda, tol, threshold, regl);
                method6 = 'LinMethod';
            end

            if MethodBool(7) == 1;
				[eps_method7, phi_method7, ampl_method7, NumberSignal7, P7, Iteration7] =...
					psi_method(Ampl_temp, Phaza_temp, X, m, lamda,eps0, limit, N0, threshold, H, h1, regl);
				method7 = 'PsiMethod';
            end

            if MethodBool(8) == 1;
                [eps_method8, phi_method8, ampl_method8, NumberSignal8, P8, Iteration8] =...
                    psi4_method(Ampl_temp, Phaza_temp, X, m, lamda,eps0,limit,N0, threshold, H, h1, regl);
                method8 = 'Psi4Method';
            end

            if MethodBool(9) == 1;
                if 1
                    [eps_method9, phi_method9 ampl_method9, NumberSignal9, P9, Iteration9] =...
                        polp_one_method(Ampl_temp, Phaza_temp, X, m, lamda, eps0, limit, N0, H, h1, regl);
                    method9 = 'PolpOne';
                else
                    [eps_method1, phi_method1, ampl_method1, NumberSignal1, P1, Iteration1] =...
                        polp_one_method(Ampl_temp, Phaza_temp, X, m, lamda, eps, limit, N0, H, h1, regl);
                    method1 = 'A';
                    [eps_method2, phi_method2, ampl_method2, NumberSignal2, P2, Iteration2] =...
                        polp_one_method(Ampl_temp, Phaza_temp, X, m, lamda, eps0, limit, N0, H, h1, regl);
                    method2 = 'B';
                    limit = [-30*pi/180 30*pi/180; -30*pi/180 30*pi/180; -30*pi/180 30*pi/180; -30*pi/180 30*pi/180];
                    [eps_method3, phi_method3, ampl_method3, NumberSignal3, P3, Iteration3] =...
                        polp_one_method(Ampl_temp, Phaza_temp, X, m, lamda, eps0, limit, N0, H, h1, regl);
                    method3 = 'C';
                end
            end

            if MethodBool(10) == 1;
				[eps_method10, phi_method10, ampl_method10, NumberSignal10, P10, Iteration10] =...
					polp_iter_method(Ampl_temp, Phaza_temp, X, m, lamda,eps0, limit, N0, threshold, H, h1, regl);
				method10 = 'PolpOneIter';
            end

            if MethodBool(11) == 1;
                [eps_method11, phi_method11, ampl_method11, NumberSignal11, P11, Iteration11] =...
                    my_method(Ampl_temp, Phaza_temp, X, m, lamda, tol, threshold, regl);
                method11 = 'MyMethod';
            end	
		
            eps_method11   = sqrt(1./(12.83*mean(Q.^2)));
            phi_method11   = sqrt(1./(12.83*mean(Q.^2)));
            ampl_method11  = sqrt(1./(12.83*mean(Q.^2)));
            NumberSignal11 = sqrt(1./(12.83*mean(Q.^2)));
            P11 = sqrt(1./(12.83*mean(Q.^2)));
            Iteration11 = sqrt(1./(12.83*mean(Q.^2)));

            %% %%%%%%%%% вычисляем ошибку измерения угла %%%%%%%%%%%%%%%%%
			% цикл по каждому методу
            for met = 1:10
                if ~strcmp('NULL',  eval(['method',int2str(met)])) % если конкретный метод был выбран

                    ogr_ = min(eval(['NumberSignal',int2str(met)]), ogr);
                    for temp = 1:length(eval(['eps_method',int2str(met)]))
                        evalin('base', [ 'Angle_', (['method',int2str(met)]),'_', int2str(temp),'(mm,zz) = eps_',...
                            (['method',int2str(met)]), '(',int2str(temp), ');']);
                    end

                    for temp = 1:length(eval(['P',int2str(met)]))
                        evalin('base', [ 'Power_', (['method',int2str(met)]),'_', int2str(temp),'(mm,zz) = P',...
                            int2str(met), '(',int2str(temp), ');']);
                    end

                    string_temp = '';
                    evalin('base', [ '[Err_', (['method',int2str(met)]), ' Sko_', (['method',int2str(met)]),...
                        '] = error_(eps, eps_', (['method',int2str(met)]), '(1:ogr_), ogr_, alpha, regl);']);
                    for temp = 1:length(eval(['P',int2str(met)]))
                        string_temp = [string_temp, ' Err_', (['method',int2str(met)]), '_', int2str(temp),...
                            '(mm,zz)', ' Sko_', (['method',int2str(met)]),  '_', int2str(temp),'(mm,zz)'];
                    end

                    string_temp = ['[', string_temp, ' ]'];
                    evalin('base', [ string_temp, ' = evalute_error_(Err_', (['method',int2str(met)]), ', Sko_',...
                        (['method',int2str(met)]),  ');']);
                    evalin('base', [ 'Iter', int2str(met), '(mm,zz) = Iteration', int2str(met), ';']);
                    evalin('base', [ 'NSignal', int2str(met), '(mm,zz) = NumberSignal', int2str(met), ';']);
                end
            end
            %%

            dT = dT + toc(Ttemp);
            T_predict = (length(Z)*Mean)*dT/( (zz-1)*Mean + mm);

            str = sprintf('Начало моделирования: %s.\n Предполагаемое дата окончания моделирования %s. \n Предполагаемое время моделирования %f сек.', datestr(T0), datestr(addtodate(datenum(T0), floor(T_predict),'second')), T_predict);
            set(Txt, 'String', str);
            waitbar( ( (zz-1)*Mean + mm) / (Mean*length(Z)))

        end
        for met = 1:12
            if ~strcmp('NULL',  eval(['method',int2str(met)]))
                for temp = 1:length(eval(['P',int2str(met)]))
                    evalin('base', [  ' SKO_method',int2str(met),  '_', int2str(temp),'(zz)=  std(Angle_method'...
                        int2str(met),  '_', int2str(temp),'(:, zz));']);
                end
            end
        end
    end


    str = sprintf(' \n Полное время моделирования %f сек.\n', dT);

    disp(str);
    close('Time')
    close(h)
    %%
	% завершение основных циклов моделирования
	

    if strcmp('SNR', method)
        Z = mean(NOISE);
    end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Выводим результаты оценок алгоритмов во вспомогательное окно

    ScreenSize = get(0, 'ScreenSize');
    otstup = 100;
	
    WindowSize = [ScreenSize(1)+2*otstup-1 ScreenSize(2)+otstup-1 ScreenSize(3)-3*otstup ScreenSize(4)-2*otstup];
    FigSize1 = [ WindowSize(3)/3-2*otstup  otstup WindowSize(3)-8*otstup WindowSize(4)-1.5*otstup];
    TextSize = [ otstup, WindowSize(4)-2.5*otstup+7, 280, 100];

    Fig1 = figure('Name', 'FAR',...
        'NumberTitle', 'Off');%,...
  
    Fig2 = figure('MenuBar', 'None',...
        'Name', 'Tune',...
        'NumberTitle', 'Off',...
        'Position', [20 200 160 330]);

    List = {'Error', 'Sko', 'NSignal', 'Angle', 'Power', 'Iteration'};
    List2 = {1, 2, 3, 'all'};

    L2 =  uicontrol('Style', 'PopupMenu','String', List2, 'Position', [10, 283, 130, 30]);
    L1 =  uicontrol('Style', 'PopupMenu','String', List, 'Position', [10, 260, 130, 30]);
    m1 = uicontrol('Style', 'RadioButton',  'String', 'MpOneMethod', 'Tag', 'MpOneMethod','Position', [10, 240, 130, 20], 'Value', MethodBool(1));
    m2 = uicontrol('Style', 'RadioButton',  'String', 'MpOneIterMethod', 'Tag', 'MpOneIterMethod','Position', [10, 220, 130, 20], 'Value', MethodBool(2));
    m3 = uicontrol('Style', 'RadioButton',  'String', 'MpMethodErm', 'Tag', 'MpMethodErm','Position', [10, 200, 130, 20], 'Value', MethodBool(3));
    m4 = uicontrol('Style', 'RadioButton',  'String', 'MpIterMethodErm', 'Tag', 'MpIterMethodErm','Position', [10, 180, 130, 20], 'Value', MethodBool(4));
    m5 = uicontrol('Style', 'RadioButton',  'String', 'LinMethodOdn', 'Tag', 'LinMethodOdn','Position', [10, 160, 130, 20], 'Value', MethodBool(5));
    m6 = uicontrol('Style', 'RadioButton',  'String', 'LinMethod', 'Tag', 'LinMethod','Position', [10, 140, 130, 20], 'Value', MethodBool(6));
    m7 = uicontrol('Style', 'RadioButton',  'String', 'PsiMethod', 'Tag', 'PsiMethod','Position', [10, 120, 130, 20], 'Value', MethodBool(7));
    m8 = uicontrol('Style', 'RadioButton',  'String', 'Psi4Method', 'Tag', 'Psi4Method','Position', [10, 100, 130, 20], 'Value', MethodBool(8));
    m9 = uicontrol('Style', 'RadioButton',  'String', 'PolpOne', 'Tag', 'PolpOne','Position', [10, 80, 130, 20], 'Value', MethodBool(9));
    m10 = uicontrol('Style', 'RadioButton',  'String', 'PolpOneIter', 'Tag', 'PolpOneIter','Position', [10, 60, 130, 20], 'Value', MethodBool(10));
    m11 = uicontrol('Style', 'RadioButton',  'String', 'Граница Крамера-Рао', 'Tag', 'MyMethod','Position', [10, 40, 130, 20], 'Value', MethodBool(11));
    p1 = uicontrol('Style', 'PushButton',  'String', 'plot', 'Tag', 'plot','Position', [10, 20, 130, 20], 'CallBack', 'output');
    p2 = uicontrol('Style', 'PushButton',  'String', 'close', 'Tag', 'close','Position', [10, 0, 130, 20], 'CallBack', 'close all');
    
	
    R = [get(m1, 'Value') get(m2, 'Value') get(m3, 'Value') get(m4, 'Value') get(m5, 'Value') get(m6, 'Value') get(m7, 'Value') get(m8, 'Value') get(m9, 'Value') get(m10, 'Value') 0 0];
    if strcmp('SNR', method)
        EdIzm = ', дБ.';
    else
        EdIzm = ', град.';
    end

    i = 0;
    A(1).sym = 'NULL';
    for tt = 1:12
        if R(tt)
            i = i +1;
            A(i).sym=eval(['method', int2str(tt)]);
        end
    end

    switch get(L1, 'Value')
        case 1
            figure(Fig1)
            clf
            set(gcf, 'Name', 'ошибка измерения луча ');
            set(gcf,'DefaultAxesLineStyleOrder', '-| --+|:o|-|--|:+',  'DefaultAxesColorOrder',[0 0 1;1 0 0;0 1 0]);
            for tt = 1:12
                if R(tt)
                    if get(L2, 'Value') < 4
                        ttt = get(L2, 'Value');
                        plot(Z, 180/pi*mean(eval(['Err_method', int2str(tt), '_',  int2str(ttt)'])));
                        hold all;
                    else
                        for ttt = 1:2
                            plot(Z, 180/pi*mean(eval(['Err_method', int2str(tt), '_',  int2str(ttt)'])));
                            hold all;
                        end
                    end
                end
            end
            xlabel([method, EdIzm], 'FontWeight', 'bold ', 'FontUnits', 'normalized');
            ylabel('Ошибка измерения угла, град.', 'FontWeight', 'bold ', 'FontUnits', 'normalized');

        case 2
            figure(Fig1)
            clf
            set(gcf, 'Name', 'ско ошибки измерения луча ');
            set(gcf,'DefaultAxesLineStyleOrder', '-| --+|:o|-|--|:+',  'DefaultAxesColorOrder',[0 0 1;1 0 0;0 1 0]);
            for tt = 1:12
                if R(tt)
                    if get(L2, 'Value') < 4
                        ttt = get(L2, 'Value');
                        plot(Z, 180/pi*((eval(['SKO_method', int2str(tt), '_',  int2str(ttt)']))));
                        hold all;
                    else
                        for ttt = 1:2
							plot(Z, 180/pi*((eval(['Sko_method', int2str(tt), '_',  int2str(ttt)']))));
                            hold all;
                        end
                    end
                end
            end
            xlabel([method, EdIzm], 'FontWeight', 'bold ', 'FontUnits', 'normalized');
            ylabel('Ско ошибки измерения угла, град.', 'FontWeight', 'bold ', 'FontUnits', 'normalized');

        case 3
            figure(Fig1)
            clf
            set(gcf, 'Name', 'колличество определенных лучей ');
            set(gcf,'DefaultAxesLineStyleOrder', '-| --+|:o|-|--|:+',  'DefaultAxesColorOrder',[0 0 1;1 0 0;0 1 0]);
            for tt = 1:12
                if R(tt)
                    plot(Z, mean(eval(['NSignal', int2str(tt)])));
                    hold all;
                end
            end
            xlabel([method, EdIzm], 'FontWeight', 'bold ', 'FontUnits', 'normalized');
            ylabel('Количество измеренных лучей', 'FontWeight', 'bold ', 'FontUnits', 'normalized');

        case 4
            figure(Fig1)
            clf
            set(gcf, 'Name', 'оценка луча ');
            set(gcf,'DefaultAxesLineStyleOrder', '-| --+|:o|-|--|:+',  'DefaultAxesColorOrder',[0 0 1;1 0 0;0 1 0]);
            for tt = 1:12
                if R(tt)
                    if get(L2, 'Value') < 4
                        ttt = get(L2, 'Value');
                        plot(Z, 180/pi*mean(eval(['Angle_method', int2str(tt), '_',  int2str(ttt)'])));
                        hold all;
                    else
                        for ttt = 1:2
                            plot(Z, 180/pi*mean(eval(['Angle_method', int2str(tt), '_',  int2str(ttt)'])));
                            hold all;
                        end
                    end
                end
            end
            xlabel([method, EdIzm], 'FontWeight', 'bold ', 'FontUnits', 'normalized');
            ylabel('Оценка угла падения волны, град', 'FontWeight', 'bold ', 'FontUnits', 'normalized');

        case 5
            figure(Fig1)
            clf
            set(gcf, 'Name', 'мощность, падающая на ФАР для каждой итерации  ');
            set(gcf,'DefaultAxesLineStyleOrder', '-| --+|:o|-|--|:+',  'DefaultAxesColorOrder',[0 0 1;1 0 0;0 1 0]);
            for tt = 1:12
                if R(tt)
                    if get(L2, 'Value') < 4
                        ttt = get(L2, 'Value');
                        plot(Z, mean(eval(['Power_method', int2str(tt), '_', int2str(ttt)])));
                        hold all;
                    else
                        for ttt = 1:2
                            plot(Z, mean(eval(['Power_method', int2str(tt), '_', int2str(ttt)])));
                            hold all;
                        end
                    end
                end
            end
            xlabel([method, EdIzm], 'FontWeight', 'bold ', 'FontUnits', 'normalized');
            ylabel('Мощность падающего поля на ФАР', 'FontWeight', 'bold ', 'FontUnits', 'normalized');


        case 6
            figure(Fig1)
            clf
            set(gcf, 'Name', 'колличество итераций ');
            set(gcf,'DefaultAxesLineStyleOrder', '-| --+|:o|-|--|:+',  'DefaultAxesColorOrder',[0 0 1;1 0 0;0 1 0]);
            for tt = 1:12
                if R(tt)
                    plot(Z, mean(eval(['Iter', int2str(tt)])));
                    hold all;
                end
            end
            xlabel([method, EdIzm], 'FontWeight', 'bold ', 'FontUnits', 'normalized');
            ylabel('колличество итераций', 'FontWeight', 'bold ', 'FontUnits', 'normalized');
    end
	
	% 11-ый метод - теоретическая граница крамера-рао
    method11 = 'Граница Крамера-Рао';
    if  get(m11, 'Value')>0
        syms lam da epsilon f1 f2 f fi1 fi2 ffi PI n pr F;
		fi = (0*da*sin(epsilon)) + (0*da*sin(pr));
        for b = 1:N-1
            fi(da, epsilon, pr) = cos(fi + (b*da*sin(epsilon)) - (b*da*sin(pr)));
        end
          fi1 = (0*da*sin(epsilon)) ; 
          fi2 =(0*da*sin(pr));
        for b = 0:N-1
            fi1(da, epsilon) = fi1 + cos(b*da*sin(epsilon));
            fi2(da, pr) = fi2 + cos(b*da*sin(pr));
        end
        fi(da, epsilon, pr) = alpha(1)*alpha(2)*fi1(da, epsilon)*fi2(da, pr);
        EN = sum(abs(Ampl_noise(:).*exp(1i*Phaza_noise(:))).^2)/2;
        f1(da, epsilon, pr) =diff(fi,2, pr);
        f2(da, epsilon, pr) =diff(f1,2, epsilon);
        f2(da, epsilon, pr) = 0.010*(f2+f1)/2;
         if strcmp('SNR', method)
            rea = f2(2*pi*d/lamda, eps(1), eps(2));
        else
            rea = f2(2*pi*d/lamda, Z*pi/180, eps(2));
        end
        rea = abs(rea);
        plot(Z, 180/pi*sqrt((1/EN./rea)./mean(Q.^2)));
        i = i + 1;
        A(i).sym='Граница Крамера-Рао';
    end
    legend(A.sym);



  
