clc;
close all;
clear;
%parametry wejściowe 
D = 0:0.01:1; %1 - 
U1 = [250, 350]; % V napięcie wejściowe 
U2 = 700; %V Napięcie wyjściowe
P2 = 9e3; %kW
R = U2^2/P2;
i1_RMS = P2./(U2*(1-D)); %RMS
i_max = P2/U2;
Delata_I = 0.08 * i1_RMS; % 8 %
i_1_MAX = (i1_RMS + Delata_I/2);
Delata_U2 = 20; % V
Tjmax = 110; % C
fs = 60e3; %częstotliowść łączeń 60 kHz
Ts = 1/fs;
i=1;
%% Topologia 1 gałęziowa
% Obliczenia dławika 

D_zakres = 0.35:0.01:0.5;
D_max = -(U1./U2)+1; 
% L = (U1(1)*D)/(Delata_I*fs);
L = max(U2*(1-D_max).*D_max)./(fs*Delata_I)
% Delta_I_obliczone = (U2*(1-D)*D)/(fs*L) % do projektu przyjęta wartość indukcyjności L(2) = 43,8 mH
% obliczenia Rezystancji wyjściowej 
% figure;
% subplot(2,1,1)
% plot(D,L,'-s')
% subplot(2,1,2)
% plot(D,Delta_I_obliczone,'-s')

C = (D(65)*Ts*U2)/(Delata_U2*R)
% %% Topologia n gałęziowa
% %parametry
% n=2;
% t_on = Ts*D;
% tau = Ts/n
% L_n = max((U2.*(1-1.*D_zakres).*1.*D_zakres)./(fs*Delata_I*1)) % l_n ustwionana wartośc maksymalną dla n = 1
% Delta_I_obliczone = (U2*(1-D).*D)./(fs*L_n)
% 
% while D(i) <= 1/2
%     n_on =1;
%     tr = (D(i)-((n_on-1)/n))*Ts;
%     d(i) = tr/tau;
%     Delta_I_2_galezi(i) = (U2*(n_on - n*D(i))*d(i))/(fs*L_n*n);
%     i = i+1;
% end
% while D(i) > 1/2
%        n_on = 2;
%        tr = (D(i)-((n_on-1)/n))*Ts;
%     d(i) = tr/tau;
%        Delta_I_2_galezi(i) = (U2*(n_on - n*D(i))*d(i))/(fs*L_n*n);
%     i=i+1;
%     if i == size(D,2)+1
%         break
%     end
% 
% end
% i = 1;
% n=4;
% while D(i) <= 1/4
%     n_on =1;
%     tr = (D(i)-((n_on-1)/n))*Ts;
%     d(i) = tr/tau;
%     Delta_I_4_galezi(i) = (U2*(n_on - n*D(i))*d(i))/(fs*L_n*n);
%     i = i+1;
% end
% while D(i) > 1/4 && D(i)<= 1/2
%        n_on = 2;
%        tr = (D(i)-((n_on-1)/n))*Ts;
%         d(i) = tr/tau;
%        Delta_I_4_galezi(i) = (U2*(n_on - n*D(i))*d(i))/(fs*L_n*n);
%     i=i+1;
% 
% 
% end
% while D(i) > 1/2 && D(i)<= 3/4
%     n_on = 3;
%     tr = (D(i)-((n_on-1)/n))*Ts;
%     d(i) = tr/tau;
%     Delta_I_4_galezi(i) = (U2*(n_on - n*D(i))*d(i))/(fs*L_n*n);
%     i=i+1;
% end
% while  D(i)> 3/4
%     n_on = 4;
%     tr = (D(i)-((n_on-1)/n))*Ts;
%     d(i) = tr/tau;
%     Delta_I_4_galezi(i) = (U2*(n_on - n*D(i))*d(i))/(fs*L_n*n);
%     i=i+1;
%     if i == size(D,2)+1
%         break
%     end
% end
% 
% %obliczenia
% 
% %wachania_i1 = Delta_I_n_galezi/i_max *100
% figure;
% hold on;
% plot(D,Delta_I_obliczone,'-x')
% plot(D,Delta_I_2_galezi,'-x')
% plot(D,Delta_I_4_galezi,'-x')
% hold off;
% legend('Przekształtnik 1 gałęziwy','Przekształtnik 2 gałęziwy','Przekształtnik 4 gałęziwy')
% xlabel('Wartość współczynnika D')
% ylabel('Wartość tętnień prądu I_1 [A]')
% axis([0.35 0.50 0 inf]);
% grid on;
%% Topologia n gałęziowa indukcyjność 
%parametry
i1_do_wykresu = P2./(U2*(1-D));
for n = 1:4
    for i = 1:n
        for j = 1:size(D,2)
            if D(1,j) >=(i-1)/n && D(1,j) <=(i)/n
                tau = Ts/n
                n_on =i;
                tr = (D(1,j)-((n_on-1)/n))*Ts;
                d(n,j) = tr/tau;
                L_n(n,j) = (U2*(n_on - n*D(1,j))*d(n,j))/(fs*Delata_I(1,j)*n);
            
            end
        end
    % disp(i)
    end
% disp(n)
end
figure;
plot(D,(L_n./i1_do_wykresu)*100);
grid on;
axis tight;
legend('Przekształtnik 1 gałęziwy','Przekształtnik 2 gałęziwy','Przekształtnik 3 gałęziwy','Przekształtnik 4 gałęziwy')
xlabel('Wartość współczynnika D')
ylabel('Wartość indukcyjności dławika L [H]')
axis([0.5 0.65 0 inf]);
L_wybrany = max(L_n(2,50:65));
%% Topologia n gałęziowa tętnienia
%parametry
i1_do_wykresu = P2./(U2*(1-D));
for n = 1:4
    for i = 1:n
        for j = 1:size(D,2)
            if D(1,j) >=(i-1)/n && D(1,j) <=(i)/n
                tau = Ts/n;
                n_on =i;
                tr = (D(1,j)-((n_on-1)/n))*Ts;
                d(n,j) = tr/tau;
                Delta_i_n(n,j) = (U2*(n_on - n*D(1,j))*d(n,j))/(fs*L_wybrany*n);
            
            end
        end
    % disp(i)
    end
% disp(n)
end
figure;
plot(D,Delta_i_n);
grid on;
axis tight;
legend('Przekształtnik 1 gałęziwy','Przekształtnik 2 gałęziwy','Przekształtnik 3 gałęziwy','Przekształtnik 4 gałęziwy')
xlabel('Wartość współczynnika D')
ylabel('Wartość tętnień prądu I_1 [A]')
axis([0.5 0.65 0 inf]);
grid on;

figure;
plot(D,(Delta_i_n./i1_do_wykresu)*100);
grid on;
axis tight;
legend('Przekształtnik 1 gałęziwy','Przekształtnik 2 gałęziwy','Przekształtnik 3 gałęziwy','Przekształtnik 4 gałęziwy')
xlabel('Wartość współczynnika D')
ylabel('wartośc tętnień prądu I_1 w [%]')
axis tight;
grid on;
I_wybrany = max(Delta_i_n(4,50:65))
%% dławik sprzężony dodatnio kształtka typ E
% k = linspace(0,1,size(D,2));
kk=0:0.05:1;
Dujemnie11=0:0.01:0.5;
Dujemnie22=0.5:0.01:1;
[k, Dujemnie1]=meshgrid(kk, Dujemnie11);
[k, Dujemnie2]=meshgrid(kk, Dujemnie22);
% for i = 1:size(k,2)
%     for n = 1:size(D,2)
%       I1_L_sprz_dodatnio(n,i) = (U2*(D(1,n)-1))*(1-2*D(1,n))/(L_wybrany*(1+k(1,i))*fs);  %dodatnie
%       I1_L_sprz_ujemne(n,i) = U2*(D(1,n)-1)*(1-2*D(1,n))/(L_wybrany*(1-k(1,i))*fs); % ujemne
%     end
% end
i1_ujemnie1=((U2*(1-2*Dujemnie1).*Dujemnie1)./(L_wybrany.*(1-k).*fs))./(P2./((1-Dujemnie1).*U2))*100;
i1_ujemnie2=((U2*(2-2*Dujemnie2).*(Dujemnie2-0.5))./(L_wybrany.*(1-k).*fs))./(P2./((1-Dujemnie2).*U2)/2)*100;

i1_dodatnie1=((U2*(1-2*Dujemnie1).*Dujemnie1)./(L_wybrany.*(1+k).*fs))./((P2./((1-Dujemnie1).*U2)))*100;
i1_dodatnie2=((U2*(2-2*Dujemnie2).*(Dujemnie2-0.5))./(L_wybrany.*(1+k).*fs))./((P2./((1-Dujemnie2).*U2)))*100;
% kształtka typ U
i1_ujemnie1n=((U2*(1-2*Dujemnie1).*Dujemnie1)./(L_wybrany.*fs))./(P2./((1-Dujemnie1).*U2))*100;
i1_ujemnie2n=((U2*(2-2*Dujemnie2).*(Dujemnie2-0.5))./(L_wybrany.*fs))./(P2./((1-Dujemnie2).*U2)/2)*100;

i1_dodatnie1n=((1-k).*(U2*(1-2*Dujemnie1).*Dujemnie1)./(L_wybrany.*(1+k).*fs))./((P2./((1-Dujemnie1).*U2)))*100;
i1_dodatnie2n=((1-k).*(U2*(2-2*Dujemnie2).*(Dujemnie2-0.5))./(L_wybrany.*(1+k).*fs))./((P2./((1-Dujemnie2).*U2)/2))*100;

% z1 = (abs(I1_L_sprz_dodatnio)./i_avg);
% z2 = (abs(I1_L_sprz_ujemne)./i_avg);
% figure;
% mesh(k,D,z1);
% axis tight;%([0 1 0 1 0 10]);
% figure;
% mesh(k,D,z2);
% axis tight;
figure(7)
mesh(Dujemnie1,k,i1_dodatnie1);
shading interp
view(3)
zlabel('{(}{\itI}{_1_p_-_p}{/\itI}{_1_a_v}{)100%}', 'FontSize', 20)
xlabel('{\itD}', 'FontSize', 20)
ylabel('{\itk}', 'FontSize', 20)
title('Sprzężenie dodatnie kształtka typu E')
%colorbar

set(gca,'fontsize',20);
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultTextFontname', 'Times New Roman')
% axis([0 1 0 1 0 8])
set(gca,'ytick',[0 0.2 0.4 0.6 0.8 1])
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1])
grid on
hold
mesh(Dujemnie2,k,i1_dodatnie2);


figure(8)
mesh(Dujemnie1,k,i1_ujemnie1);
shading interp
view(3)
zlabel('{(}{\itI}{_1_p_-_p}{/\itI}{_1_a_v}{)100%}', 'FontSize', 20)
xlabel('{\itD}', 'FontSize', 20)
ylabel('{\itk}', 'FontSize', 20)
axis([0 1 0 1 0 20])
title('Sprzężenie ujemne kształtka typu E')
%colorbar

set(gca,'fontsize',20);
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultTextFontname', 'Times New Roman')
% axis([0 1 0 1 0 8])
set(gca,'ytick',[0 0.2 0.4 0.6 0.8 1])
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1])
grid on
hold
mesh(Dujemnie2,k,i1_ujemnie2);

figure(9)
mesh(Dujemnie1,k,i1_dodatnie1n);
shading interp
view(3)
zlabel('{(}{\itI}{_1_p_-_p}{/\itI}{_1_a_v}{)100%}', 'FontSize', 20)
xlabel('{\itD}', 'FontSize', 20)
ylabel('{\itk}', 'FontSize', 20)
title('Sprzężenie dodatnie kształtka typu U')
%colorbar

set(gca,'fontsize',20);
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultTextFontname', 'Times New Roman')
% axis([0 1 0 1 0 8])
set(gca,'ytick',[0 0.2 0.4 0.6 0.8 1])
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1])
grid on
hold
mesh(Dujemnie2,k,i1_dodatnie2n);

figure(10)
mesh(Dujemnie1,k,i1_ujemnie1n);
shading interp
view(3)
zlabel('{(}{\itI}{_1_p_-_p}{/\itI}{_1_a_v}{)100%}', 'FontSize', 20)
xlabel('{\itD}', 'FontSize', 20)
ylabel('{\itk}', 'FontSize', 20)
%colorbar
title('Sprzężenie ujemne kształtka typu U')

set(gca,'fontsize',20);
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultTextFontname', 'Times New Roman')
% axis([0 1 0 1 0 20])
axis([0 1 0 1 0 20])
set(gca,'ytick',[0 0.2 0.4 0.6 0.8 1])
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1])
grid on
hold
mesh(Dujemnie2,k,i1_ujemnie2n);

% %% dławik sprzężony dodatnio kształtka typ U
% k = linspace(0,1,size(D,2));
% for i = 1:size(k,2)
%     for n = 1:size(D,2)
%       I1_L_sprz_dodatnio_U(n,i) = U2*D(1,n)*(1-2*D(1,n))/(L_wybrany*fs);  %dodatnie
%       I1_L_sprz_ujemne_U(n,i) = U2*D(1,n)*(1-k(1,i))*(1-2*D(1,n))/(L_wybrany*(1+k(1,i))*fs); % ujemne
%     end
% end
% 
% z1_U = (abs(I1_L_sprz_dodatnio_U)./i_avg);
% z2_U = (abs(I1_L_sprz_ujemne_U)./i_avg);
% figure;
% mesh(D,k,z1_U);
% axis tight;%([0 1 0 1 0 10]);
% figure;
% mesh(D,k,z2_U);
% axis tight;


%% projektowanie dławika 
J = 500;
k = 0.5;
Bmax = 0.3;
i_l_MAX = 13.86; % A
I_l_RMS = 9.37; % A
Ap_e = (L_wybrany*i_l_MAX*I_l_RMS*10^4)/(Bmax*J*k);
% dobór rdzenia
A_c = [692, 645, 392, 560, 211, 447,411, 676, 840, 10.1, 100, 234,142]; % 1 rdzeń E100/60/28 , U100/57/25 + I100/25/25, E80/38/20, U126/91/20,EC70/34/17,U93/76/16,E55/28/25,E71/33/32, UI93/76/30 ,E13/6/3 , E35/18/10,E42/33/20 ,E41/17/12
W_a = [(73.15-27.5)/2*46.85*2, 25.4*25.4, (59.1-19.8)/2*28.2*2,  70*63, (44.5-16.4)/2*22.75*2, 27.5*30, (37.5-17.2)/2*18.5*2, (48-22)/2*21.9*2, 34.6*48, (9.5-3.2)/2 * 4.1*2, (24.5-10)/2 * 12.5*2,(29.5-12.2)/2 *26*2 , (28.6-12.45)/2 *10.4*2];
Ap_p = A_c.*W_a/10^4;
ilosc = Ap_p./Ap_e
% ilość zwojów
S_fe= [ (71.7-44.5)/2*16.4, (42-29.5)/2 * 20, 2*(40.6-28.6)/2 *12.4 ]% EC70/34/17, E42/33/20, E41/17/12
N = (L_wybrany*i_l_MAX)./(Bmax*S_fe/10^6)
S_cu = I_l_RMS/5; % pole przekroju poprzecznego przeewodu elektrycznego
d_cu = 2*sqrt(S_cu/pi) % mm
glebokosc_wnikania = sqrt((1.72*10^-8)/(pi*1*4*pi*10^-7*fs)); % w 0.x mm
maksymalna_srednica_przewowdu = 2*glebokosc_wnikania
pole_przekroju_licy = 2*2;
ilosc_licy = pole_przekroju_licy*N
czy_sie_zmiesci = W_a(1,[5 12 13])./ilosc_licy % wybrany rdzeń 12 E42/33/20
% wybrany rdzeń E71/33/32(8) => U93/76/16
%% sprawdzenie elementów półprzewodnikowych 
% tranzystor 
n = 4;
delta_l_4 = 9.864;
I2_av = i_max;
IT_av = D(1,51:65).*i1_RMS(1,51:65)./n;
IT_max = i1_RMS(1,51:65)/n+delta_l_4/2;
IT_rms = sqrt(D(1,51:65).*(IT_max.^2+delta_l_4^2./3-IT_max*delta_l_4));
% dioda
ID_av = (1-D(1,51:65)).*i1_RMS(1,51:65)./n;
ID_max = i1_RMS(1,51:65)/n+delta_l_4/2;
ID_rms = sqrt((1-D(1,51:65)).*(ID_max.^2+delta_l_4.^2./3-ID_max.*delta_l_4));
%współczynnik zapassu 
k = 1.5;
% szykane parametry tranzystora i diody
disp('średni prąd tranzystora z zapasem 1.5');
disp(max(IT_av)*1.5);
disp('skuteczny tranzystora prąd z zapasem 1.5');
disp(max(IT_rms)*1.5);
disp('maksymlany tranzystora prąd z zapasem 1.5');
disp(max(IT_max)*1.5);
disp('średni prąd diody z zapasem 1.5');
disp(max(ID_av)*1.5);
disp('skuteczny diody prąd z zapasem 1.5');
disp(max(ID_rms)*1.5);
disp('maksymlany diody prąd z zapasem 1.5');
disp(max(ID_max)*1.5);
IT_av = max(IT_av)
IT_max = max(IT_max)
IT_rms = max(IT_rms)
% dioda
ID_av = max(ID_av)
ID_max = max(ID_max)
ID_rms = max(ID_rms)
%% straty mocy
% straty przewodzenia
% parametry tr1 mosfet
Tr1_rdson = 1.15*120e-3;

Tr1_Pc_t = Tr1_rdson*IT_rms^2 % W
I_TR1 = 15;
U_TR1 = 700;
RG = 2.5+16;
Ciss = 350e-12;
Crss = 3e-12;
Ugg = 15;
Ugsth = 2.1;
Ugsp = 3.5;
Uds = 15;
U0 = 700;
% Coss = 1060e-12;
Eoss = 9e-6;
tri = RG * Ciss * log((Ugg - Ugsth) / (Ugg - Ugsp));
tfu = RG * Crss * (U0 / (Ugg - Ugsp));
tfi = RG * Ciss * log(Ugsp / Ugsth);
tru = RG * Crss * (Uds / Ugsp);
Eon = 0.5*U2*IT_rms*(tri+tfu)+ Eoss;
Eoff = 0.5*U2*IT_rms*(tru+tfi);
K_on = Eon/22
K_off = Eoff/22


Tr1_P_swt = fs*(Eon+Eoff)*(IT_rms/I_TR1)*(U2/U_TR1) %W
Tr1_P_loss = Tr1_P_swt+Tr1_Pc_t
% parametry tr2 mosfet
Tr2_rdson = 2.2*350e-3;

Tr2_Pc_t = Tr2_rdson*IT_rms^2; % W
RG = 10;
Ciss = 10.9e-9;
Crss = 67e-12;
Ugg = 10;
Ugsth = 3.5;
Ugsp = 6.5;
Uds = 0.5*1e3;
U0 = 1000;
Coss = 1060e-12;
Eoss = Coss*U2^2/2
tri = RG * Ciss * log((Ugg - Ugsth) / (Ugg - Ugsp));
tfu = RG * Crss * (U0 / (Ugg - Ugsp));
tfi = RG * Ciss * log(Ugsp / Ugsth);
tru = RG * Crss * (Uds / Ugsp);
Eon = 0.5*U2*IT_rms*(tri+tfu)+ Eoss;
Eoff = 0.5*U2*IT_rms*(tru+tfi);
I_TR2 = 15;
U_TR2 = 700;

Tr2_P_swt = fs*(Eon+Eoff)*(IT_rms/I_TR2)*(U2/U_TR2) %W
Tr2_P_loss = Tr2_P_swt+Tr2_Pc_t
%tr3 igbt
Uce0 = 3.2;
Rt = 50e-3;
Tr3_Pc_t = Uce0*IT_av*Rt*IT_rms^2
Eon = 2800e-6;
Eoff = 750e-6;
I_TR3 = 20;
U_TR3 = 800;

Tr3_P_swt = fs*(Eon+Eoff)*(IT_rms/I_TR3)*(U2/U_TR3) %W
Tr3_P_loss = Tr3_P_swt+Tr3_Pc_t
%dioda
Uf0 = 1.4;
Rd = 0.022;
Pcd = Uf0*ID_av*Rd*ID_rms^2
Erec = 564e-9;
I_D = 8;
U_D = 600;
% Pswd = fs*(Erec)*(ID_rms/I_D)*(U2/U_D)
If = ID_av;
trr = 47e-9/0.9;
eoffd = (U2)*0.5*If*trr/8
Pswd = fs*(U2)*0.5*If*trr
D1_P_loss = Pswd+Pcd
% dioda 2

Uf0 = 1.9;
Rd = 0.09;
Pcd = Uf0*ID_av*Rd*ID_rms^2
Ec = 30.51e-6;

Pswd = fs*Ec
D2_P_loss = Pswd+Pcd
%% estymacja radiatora
CSPI = 3
Rth = 0.15:0.05:0.8;
Vth = 1./(CSPI*Rth)