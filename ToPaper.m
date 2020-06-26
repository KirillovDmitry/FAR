% скрипт, предназначенный для построения графики для статьи
clear all; close all; clc;
setup;

[Ampl_clear, Phaza_clear, Ampl_noise, Phaza_noise, Phaza_lin, Pow, sigma2] = field(eps, alpha, phi0, lamda, M, SNR, X);

S = Ampl_clear.*exp(1i.*Phaza_clear) - alpha(1)*exp(1i*2*pi*sin(eps(1))/lamda* X);
Ampl = abs( S );
Phaza2 = angle( S );

R1 = (2*pi)/lamda*sin(eps(1))*X;
R2 = (2*pi)/lamda*sin(eps(2))*X;


p = plot(X, unwrap(Phaza_clear)/1, 'r', X, unwrap(Phaza2)/1, 'b');
set(p,'LineWidth', 3)
hold on
t = plot( X, R1/1, '--', X, R2/1, '--');
set(t, 'Color', 'black')
set(gca,'YTick',-3*pi:pi:3*pi)
set(gca,'YTickLabel',{'-3pi', '-2pi','-pi','0','pi','2pi','3pi'})
