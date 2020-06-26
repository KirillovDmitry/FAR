% скрипт для обработки и вывода результата в окно
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
            %   legend(A.sym);
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
                        %plot(Z, 180/pi*sqrt(mean(eval(['Sko_method', int2str(tt), '_',  int2str(ttt)']))));
                        plot(Z, 180/pi*((eval(['SKO_method', int2str(tt), '_',  int2str(ttt)']))));
                        hold all;
                    else
                        for ttt = 1:2
                            %plot(Z, 180/pi*sqrt(mean(eval(['Sko_method', int2str(tt), '_',  int2str(ttt)']))));
                            plot(Z, 180/pi*((eval(['Sko_method', int2str(tt), '_',  int2str(ttt)']))));
                            hold all;
                        end
                    end
                end
            end
            %legend(A.sym);
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
            %       legend(A.sym);
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
            %     legend(A.sym);
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
            %  legend(A.sym);
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
            % legend(A.sym);
            xlabel([method, EdIzm], 'FontWeight', 'bold ', 'FontUnits', 'normalized');
            ylabel('колличество итераций', 'FontWeight', 'bold ', 'FontUnits', 'normalized');

    end

    method11 = 'Граница Крамера-Рао';
    if  get(m11, 'Value')>0
        %         syms lam da epsilon f2 f PI n;
        %         f(lam, da, PI, epsilon) = cos(-2*PI*da/lam*sin(epsilon));
        %         f(da, epsilon) = cos(da*sin(epsilon)) ;
        %         hold all;
        %         f2=diff(f,2, epsilon);
        %         V = 0;
        %         for b = 0 :N-1
        %             V = V + 2*pi*b*d/lamda;
        %         end
        %         plot(Z, 180/pi*sqrt(1/V^2./mean(Q.^2)));
        %
        syms lam da epsilon f1 f2 f fi1 fi2 ffi PI n pr F;
        %f(lam, da, PI, epsilon) = cos(-2*PI*da/lam*sin(epsilon));
        %         f(da, epsilon) = cos(0*sin(epsilon))+cos(da*sin(epsilon))+cos(2*da*sin(epsilon))+cos(3*da*sin(epsilon))+cos(4*da*sin(epsilon))+cos(5*da*sin(epsilon)) -...
        %             cos(0*sin(0))+cos(da*sin(0))+cos(2*da*sin(0))+cos(3*da*sin(0))+cos(4*da*sin(0))+cos(5*da*sin(0));
        %cos( (0*da*sin(epsilon))+ (da*sin(epsilon))+ (2*da*sin(epsilon))+ (3*da*sin(epsilon))+ (4*da*sin(epsilon))+ (5*da*sin(epsilon)) - (0*sin(0)) - (da*sin(0)) - (2*da*sin(0)) - (3*da*sin(0)) - (4*da*sin(0)) - (5*da*sin(0)));
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
        % fi(da, epsilon, pr) = cos(fi(da, epsilon, pr));
        %f(da, epsilon) = cos( da*sin(epsilon));
       % f2(da, epsilon, pr) =diff(fi,2, epsilon);
        %plot(Z, 180/pi*sqrt(-EN/16/f2(2*pi*d/lamda,eps(1),eps(2))./mean(Q.^2)));
        %          plot(Z, 180/pi*sqrt(-EN/8/f2(2*pi*d/lamda,0,0)./mean(Q.^2)));
        %          Z(zz)*pi/180;
        % WWW = sqrt(-(EN/8)./double(f2(2*pi*d/lamda, Z*pi/180, Z*pi/180))./mean(Q.^2));
        %          plot(Z, 180/pi*(real(WWW) + 0.5*imag(WWW)));
        f1(da, epsilon, pr) =diff(fi,2, pr);
        f2(da, epsilon, pr) =diff(f1,2, epsilon);
        f2(da, epsilon, pr) = 0.010*(f2+f1)/2;
        %f2(da, epsilon, pr) = 0.00010*(f2+f1)/2;
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



