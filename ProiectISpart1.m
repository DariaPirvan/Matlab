t=Pirvan(:,1)
u=Pirvan(:,2)
y=Pirvan(:,3)
plot(t,u,t,y) %A100 D1003
%ce este semnalul de intrare? - sin - amplitudinea 1
%defazaj - intarziereea dintre 2semnale in care cele 2 semnale trec prinn
%acc val - 
%u isi schimba frecventa
%%
%ce este rezonanta? amplitudini dif ale semnalelor de iesire
% de scris care e eroarea
i1=270;
i2=293;
u0=mean(u(i1:i2))
y0=mean(y(i1:i2))
i3=266;
i4=287;
Mr = (y(i1)-y(i2))/(u(i3)-u(i4))%amplificare maxima a sistemului( la rezonanta) -> ampl semnalului de iesire/ampl semnalului de intrare
Tr= 2* (t(i2)-t(i1))%perioada semnalului de iesire la rezonanta -: pt a calc wr
wr = (2*pi)/Tr %pulsatia de rezonanta
zetta = sqrt((Mr-sqrt(Mr^2-1))/2/Mr)%zetta = sqrt(Mr - sqrt(Mr^2-1))/(2*Mr)%zetta=sqrt((Mr-sqrt(Mr^2-1))/2*Mr) %trebuie sa fie mai mic decat rad2/2
wn = wr/sqrt(1-2*zetta^2) %legatura intre wr si wn este prin zetta!
num = wn^2
den = [1 2*zetta*wn wn^2]
% k=(y(i3)-y(i4))/(u(i1)-u(i2))
k=mean(y)/mean(u)
%k=1; %il pot afla si de la inceputul semnalului?!!
%spatiul starilor
A=[0,1;-wn^2,-2*zetta*wn]
B=[0;k*wn^2]
C=[1,0]
D=[0]
sys=ss(A,B,C,D);
ysim=lsim(sys,u,t,[y(1),(y(2)-y(1))/(t(i2)-t(i1))]);
figure
plot(t,y,t,ysim)
J=norm(y-ysim)/sqrt(length(y))
%eroarea medie patratica?
eMPN=norm(y-ysim)/norm(y-mean(y))

%% identific pentru semnalul de intrare
Tr= 2* (t(i4)-t(i3))
wr = (2*pi)/Tr
wn = wr/ sqrt(1-2*zetta^2)
k=(y(i3)-y(i4))/(u(i1)-u(i2)) %??
k=mean(y)/mean(u)
A=[0,1;-wn^2,-2*zetta*wn]
B=[0;k*wn^2]
C=[1,0]
D=[0]
sys=ss(A,B,C,D);
ysim=lsim(sys,u,t,[y(1),(y(2)-y(1))/(t(i2)-t(i1))]);
figure
plot(t,y,t,ysim)
J=norm(y-ysim)/sqrt(length(y))
%eroarea medie patratica?
eMPN=norm(y-ysim)/norm(y-mean(y))

% num = wn^2
% den = [1 2*zetta*wn wn^2]
% yc = lsim(num, den, u, t) 
% plot(t,[y,yc]) 


%% Identificarea Nyquist
%wn pe nyquist la intersectia cu Oy
%wn pe bode e la defazaj 90 grade
%defazajul e scaderea dintre doi X de pe sistemul initial
%in ce se masoara defazajul? - nu in radiani
nyquist(num,den)
dt = t(i2)-t(i4) %i1 i3
%cum calc defazajul de  la ieseire? in fctie de ce
ph1= dt* wr %defazajul fata de ce? in radiani
ph1 = (dt* wr * 180)/pi %transform din radiani in grade
i5 = 434;
i6 = 450;
i7= 425;
i8= 442;
dt = t(i5)-t(i7) %iesire max - intrare max
TN = 2* (t(i8)-t(i7)) %iesire max - iesire min sau intrare
wN1 = (2*pi)/TN
ph = (dt* wN1 * 180)/pi %ph = 95.2941?
M1 = (y(i5)-y(i6))/(u(i7)-u(i8)) %amplitud intr / ampl iesire
zetta1 = 1/M1/2
num1 = wN1^2
den1 = [1 2*zetta1*wN1 wN1^2]
% yc1 = lsim(num1, den1, u, t) ;
% figure
% plot(t,[y,yc1])
% k=1;
k=(y(i7)-y(i8))/(u(i5)-u(i6))
k=mean(y)/mean(u)
A=[0,1;-wN1^2,-2*zetta1*wN1]
B=[0;k*wN1^2]
C=[1,0]
D=[0]
sys=ss(A,B,C,D);
ysim=lsim(sys,u,t,[y(1),(y(2)-y(1))/(t(i2)-t(i1))]);
figure
plot(t,y,t,ysim)
figure
nyquist(num1, den1)
title('Nyquist 2')
J=norm(y-ysim)/sqrt(length(y))
%eroarea medie patratica?
eMPN=norm(y-ysim)/norm(y-mean(y))


%% cu D Bode
figure
bode(num1,den1)
wr1 = wN1 * sqrt(1-2*zetta1^2)
i9=854;  %pe albastru 854 865 848 858
i10=862;
i11=848;
i12=856;
T2 = 2* (t(i10)-t(i9))

%pe bode la faza de -90, luam frecventa si o cautam in modul si luam de pe grafic
% modulul in dB
k=1
modul=k/(2*zetta1) %il verif in Nyquist
modul_abs=10^(1.2759/20) %val absoluta a modulului 

M2=(y(i9)-y(i10))/(u(i11)-u(i12)) 
M2db = 20*log10(M2)%??

