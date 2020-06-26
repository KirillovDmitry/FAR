% функция нахождения минимума функционала ошибки измерения угловой координаты по расскрыву ФАР
function [xval, yval, output] = my_fminsearch(Func,ampl, phi, X, lamda, Y, limit, N0, H, h)
    IT = 0;
    S = size(limit);

    if S(1)==1
        Yval = zeros(S(1), N0);
        Xval = zeros(S(1), N0);
        range1 = linspace(limit(1,1), limit(1, 2), N0+1);
        for i = 1:N0
            x0 = (range1(i+1) +  range1(i))/2;
            switch (func2str(Func))
                case 'mp_one'
                    [Xval(i), Yval(i),  exitflag, out]  =  fminsearch(@(xx)mp_one(xx, ampl(1), X, lamda, Y),  x0,  optimset(H, h(1), 'MaxIter', h(2)));
                    msgstr = lastwarn;
                    if (Yval(i)<-100) || (exitflag~=1) || ( ~isempty(msgstr)) || ~( (limit(1,1)<Xval(i) && Xval(i)<limit(1,2)) )
                        Yval(i) = 1;
                        Xval(i) = 0;
                        lastwarn('')
                    end
                    IT = IT + out.iterations;
                case 'polp'
                    [Xval(i), Yval(i),  exitflag, out]  =  fminsearch(@(xx)polp(xx, ampl(1), X, lamda, Y),  x0,  optimset(H, h(1), 'MaxIter', h(2)));
                    msgstr = lastwarn;
                    if (exitflag~=1) || ( ~isempty(msgstr)) || ~( (limit(1,1)<Xval(i) && Xval(i)<limit(1,2)) )
                        Yval(i) = 1;
                        Xval(i) = 0;
                        lastwarn('')
                    end
                    IT = IT + out.iterations;
                case 'kepon_one'
                    [Xval(i), Yval(i),  exitflag, out]  =  fminsearch(@(xx)polp(xx, ampl(1), X, lamda, Y),  x0,  optimset(H, h(1), 'MaxIter', h(2)));
                    msgstr = lastwarn;
                    if (Yval(i)<-100) || (exitflag~=1) || ( ~isempty(msgstr)) || ~( (limit(1,1)<Xval(i) && Xval(i)<limit(1,2)) )
                        Yval(i) = 1;
                        Xval(i) = 0;
                        lastwarn('')
                    end
                    IT = IT + out.iterations;
            end
        end

        [C, I] = min(Yval);
        yval = C;
        if yval > 0
            xval = 0;
        else
            xval = Xval(I);
        end

    else
        range1 = linspace(limit(1,1), limit(1, 2), N0+1);
        range2 = linspace(limit(2,1), limit(2, 2), N0+1);

        Yval = zeros(N0, N0);
        X1 = zeros(N0, N0);
        X2 = zeros(N0, N0);
        X3 = zeros(N0, N0);
        X4 = zeros(N0, N0);
        for i = 1:N0
            for j = 1:N0
                x1 = (range1(i+1) +  range1(i))/2;
                x2 = (range2(j+1) +  range2(j))/2;
                switch (func2str(Func))
                    case 'mp_erm'
                        [temp_ Yval(i,j) exitflag out]  =  fminsearch(@(xx)mp_erm(xx, ampl(1:2), X, lamda, Y),  [x1 x2],  optimset(H, h(1), 'MaxIter', h(2)));
                        msgstr = lastwarn;
                        if (Yval(i,j)<-100) || (exitflag~=1) || ( ~isempty(msgstr)) || ~( (limit(1,1)<temp_(1) && temp_(1)<limit(1,2)) ) || ~( (limit(2,1)<temp_(2) && temp_(2)<limit(2,2)) )
                            Yval(i,j) = 1;
                            X1(i,j) = 0;
                            X2(i,j) = 0;
                            lastwarn('')
                        else
                            X1(i,j) = temp_(1);
                            X2(i,j) = temp_(2);
                        end
                        IT = IT + out.iterations;
                    case 'psi'
                        [temp_ Yval(i,j) exitflag out]  =  fminsearch(@(xx)psi(xx, phi, ampl(1:2), X, lamda, Y),  [x1 x2],  optimset(H, h(1), 'MaxIter', h(2)));
                        msgstr = lastwarn;
                        if (Yval(i,j)>100) || (exitflag~=1) || ( ~isempty(msgstr)) || ~( (limit(1,1)<temp_(1) && temp_(1)<limit(1,2)) ) || ~( (limit(2,1)<temp_(2) && temp_(2)<limit(2,2)) )
                            Yval(i,j) = 10000;
                            X1(i,j) = 0;
                            X2(i,j) = 0;
                            lastwarn('')
                        else
                            X1(i,j) = temp_(1);
                            X2(i,j) = temp_(2);
                            Yval(i,j) = - 1/Yval(i,j);
                        end
                        IT = IT + out.iterations;
                    case 'psi4'
                        [temp_ Yval(i,j) exitflag out]  =  fminsearch(@(xx)psi4(xx, phi,X, lamda, Y),  [x1 x2 ampl(1:2)],  optimset(H, h(1), 'MaxIter', h(2)));
                        msgstr = lastwarn;
                        if (Yval(i,j)>100) || (exitflag~=1) || ( ~isempty(msgstr)) || ~( (limit(1,1)<temp_(1) && temp_(1)<limit(1,2)) ) || ~( (limit(2,1)<temp_(2) && temp_(2)<limit(2,2)) )
                            Yval(i,j) = 10000;
                            X1(i,j) = 0;
                            X2(i,j) = 0;
                            X3(i,j) = 0;
                            X4(i,j) = 0;
                            lastwarn('')
                        else
                            X1(i,j) = temp_(1);
                            X2(i,j) = temp_(2);
                            X3(i,j) = temp_(3);
                            X4(i,j) = temp_(4);
                            Yval(i,j) = - 1/Yval(i,j);
                        end
                        IT = IT + out.iterations;
                    case 'psi45'
                        [temp_ Yval(i,j) exitflag out]  =  fminsearch(@(xx)psi45(xx, X, lamda, Y),  [x1 x2 ampl(1) phi(1)],  optimset(H, h(1), 'MaxIter', h(2)));
                        msgstr = lastwarn;
                        if (Yval(i,j)>100) || (exitflag~=1) || ( ~isempty(msgstr)) || ~( (limit(1,1)<temp_(1) && temp_(1)<limit(1,2)) ) || ~( (limit(2,1)<temp_(2) && temp_(2)<limit(2,2)) )
                            Yval(i,j) = 10000;
                            X1(i,j) = 0;
                            X2(i,j) = 0;
                            X3(i,j) = 0;
                            X4(i,j) = 0;
                            lastwarn('')
                        else
                            X1(i,j) = temp_(1);
                            X2(i,j) = temp_(2);
                            X3(i,j) = temp_(3);
                            X4(i,j) = temp_(4);
                            Yval(i,j) = - 1/Yval(i,j);
                        end
                        IT = IT + out.iterations;
                    case 'psi6'
                        [temp_ Yval(i,j) exitflag out]  =  fminsearch(@(xx)psi4(xx, X, lamda, Y),  [x1 x2 ampl(1:2) phi(1:2)],  optimset(H, h(1), 'MaxIter', h(2)));
                        msgstr = lastwarn;
                        if (Yval(i,j)>100) || (exitflag~=1) || ( ~isempty(msgstr)) || ~( (limit(1,1)<temp_(1) && temp_(1)<limit(1,2)) ) || ~( (limit(2,1)<temp_(2) && temp_(2)<limit(2,2)) )
                            Yval(i,j) = 10000;
                            X1(i,j) = 0;
                            X2(i,j) = 0;
                            X3(i,j) = 0;
                            X4(i,j) = 0;
                            X5(i,j) = 0;
                            X6(i,j) = 0;
                            lastwarn('')
                        else
                            X1(i,j) = temp_(1);
                            X2(i,j) = temp_(2);
                            X3(i,j) = temp_(3);
                            X4(i,j) = temp_(4);
                            X5(i,j) = temp_(5);
                            X6(i,j) = temp_(6);
                            Yval(i,j) = - 1/Yval(i,j);
                        end
                        IT = IT + out.iterations;
                end

            end
        end

        [~, I1] = min(Yval);
        [C2, I2] = min(min(Yval));
        yval = C2;
        if yval > 0
            switch (func2str(Func))
                case 'mp_erm'
                    xval = 0;%fminsearch(@(xx)mp_erm(xx, ampl(1:2), X, lamda, Y),  [1*180/pi -1*180/pi],  optimset(H, h(1), 'MaxIter', h(2)));
                case 'psi'
                    xval = 0;%fminsearch(@(xx)psi(xx, phi, ampl(1:2), X, lamda, Y),  [1*180/pi -1*180/pi],  optimset(H, h(1), 'MaxIter', h(2)));
                case 'psi4'
                    xval = 0;%fminsearch(@(xx)psi4(xx, phi, X, lamda, Y),  [1*180/pi -1*180/pi 1 0.9],  optimset(H, h(1), 'MaxIter', h(2)));
                case 'psi45'
                    xval = 0;
                case 'psi6'
                    xval = 0;
            end
        else
            xval(1) = X1(I1(I2),I2);
            xval(2) = X2(I1(I2),I2);
            switch (func2str(Func))
                case 'psi4'
                    xval(3) = X3(I1(I2),I2);
                    xval(4) = X4(I1(I2),I2);
                case 'psi45'
                    xval(3) = X3(I1(I2),I2);
                    xval(4) = X4(I1(I2),I2);
                case 'psi6'
                    xval(3) = X3(I1(I2),I2);
                    xval(4) = X4(I1(I2),I2);
                    xval(5) = X5(I1(I2),I2);
                    xval(6) = X6(I1(I2),I2);
            end
        end
    end

    output.iterations = IT;
