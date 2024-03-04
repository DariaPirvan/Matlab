t=Pirvan(:,1) %variatia in timp a unei intrari si a doua iesiri
u=Pirvan(:,2)  %intrarea
y1=Pirvan(:,3) %prima iesire
y2=Pirvan(:,4) %a doua iesire
plot(t,u,t,y1), title('Pentru prima iesire')
figure
plot(t,u,t,y2), title('Pentru a doua iesire')
%%
%de pe semnalul de iesire
i1=368;%452;%270;
i2=391;%473;%293;
%de pe semnalul de intrare
i3=363;%446;%266;
i4=384%468;%287;

k=mean(y1)/mean(u)%factor de proportionalitate
Mr = (y1(i1)-y1(i2))/(u(i3)-u(i4))/k %amplificare maxima a iesirii la rez
Tr= 2* (t(i4)-t(i3))%perioada la rezonanta a semnalului de intrare                   -: pt a calc wr
wr = (2*pi)/Tr%pulsatia de rezonanta
zetta = sqrt((Mr-sqrt(Mr^2-1))/2/Mr)                                                     %trebuie sa fie mai mic decat rad2/2
wn = wr/sqrt(1-2*zetta^2)%legatura intre wr si wn este prin zetta!
% num = wn^2
% den = [1 2*zetta*wn wn^2]
% k=(y(i3)-y(i4))/(u(i1)-u(i2))
%spatiul starilor
A=[0,1;-wn^2,-2*zetta*wn]
B=[0;k*wn^2]
C=[1,0]
D=[0]
ysim=lsim(A,B,C,D,u,t,[y1(1),(y1(2)-y1(1))/(t(2)-t(1))]); %dt- perioada de esantionare      %dt = t(i2)-t(i1) 
figure
plot(t,y1,t,ysim)
%eroarea medie patratica
J=norm(y1-ysim)/sqrt(length(y1))
%eroarea medie patratica normalizata
eMPN=norm(y1-ysim)/norm(y1-mean(y1))

%% Identificarea Nyquist
%wn pe nyquist la intersectia cu Oy
%wn pe bode e la defazaj 90 grade
% figure
% nyquist(num,den), title("lslsl")
% dt = t(i2)-t(i4)
% ph1= dt* wr 
% ph1 = (dt* wr * 180)/pi
i5 = 469;
i6 = 486;
i7= 460;
i8= 477;

dt = t(i7)-t(i5) %iesire max - intrare max %INTARZIEREA
TN = 2* (t(i8)-t(i7)) %iesire max - iesire min sau intrare
wN1 = (2*pi)/TN
ph = (dt* wN1 * 180)/pi %ph = 95.2941?
M1 = (y1(i5)-y1(i6))/(u(i7)-u(i8))/k %amplitud intr / ampl iesire
zetta1 = k/M1/2
num1 = k*wN1^2
den1 = [1 2*zetta1*wN1 wN1^2]
k=mean(y1)/mean(u);%sau din FCJ de pe Nyquist
A=[0,1;-wN1^2,-2*zetta1*wN1]
B=[0;k*wN1^2]
C=[1,0]
D=[0]
sys=ss(A,B,C,D);
ysim=lsim(sys,u,t,[y1(1),(y1(2)-y1(1))/(t(i2)-t(i1))]);
figure
plot(t,y1,t,ysim)
figure
nyquist(num1, den1) %frecv trebe sa fie egala cu cea calc mai sus si img trebe sa fie egala cu modul de M1
title('Nyquist')
J=norm(y1-ysim)/sqrt(length(y1))
eMPN=norm(y1-ysim)/norm(y1-mean(y1))

%% verificare cu D Bode
figure
bode(num1,den1)
wr1 = wN1 * sqrt(1-2*zetta1^2) %wr1 pctul max de pe modul bode
%indici de la sfarsitul semnalului
i9=854;
i10=862;
i11=848;
i12=856;
T2 = 2* (t(i10)-t(i9)) %trebuie ca T2 =TN/2 Tn = 0.0017 si T2 = 0.0008 - da bine

%pe bode la faza de -90, luam frecventa si o cautam in modul si luam de pe grafic
% modulul in dB
% k=1
% modul=k/(2*zetta1) %il verif in Nyquist
% modul_abs=10^(1.2759/20) %val absoluta a modulului 

M2=(y1(i9)-y1(i10))/(u(i11)-u(i12))/k
M2db = 20*log10(M2)% amplificarea in db -> trebuie sa dea -12(pe octava)->am panta de -40

%% Identificare parametrica
dt = t(2)-t(1) %pas de achizitie
d_id = iddata(y1,u,dt)
%model calculat cu arx (identificare prin arx)
Marx = arx(d_id,[2,2,1]) %model calculat cu arx
Hz = tf(Marx.B, Marx.A,dt) %fctia de transfer in discret
Hs = d2c(Hz,'zoh') % fctia de tr ansfer in continuu %apare zero de la 'zoh'
% Hs = tf(1.534e07,[1 4610 1.505e07]);
figure;resid(d_id, Marx,5) %stg autocorelatie/ dr intercorelatie - 5 corelatii
figure;compare(d_id, Marx) 

%la arx nu trece testul(pt ca nu trece la autocorelatie, intercorelatia
%trece dar nu ne intereseaza)
%% Rafinare prin post procesare(pem) pt arx
%iau cea mai rea identificare(arx in cazul meu) si fac cu metoda err de predictie(pem)
Marx_pem=pem(Marx,d_id)
%validare statistica
resid(d_id,Marx_pem)
figure
%gradul de suprapunere
compare(d_id,Marx_pem)
%% Identif cu variab instrumentale %nu iese intercorelatia, si nici autocor (dar nu ne intereseaza)
Mvi = iv4(d_id, [2,2,1])
figure;resid(d_id, Mvi, 5)
figure;compare(d_id,Mvi) %aici am comparat Marx cu Mvi, Mvi - 71.162% mai prost ca arx

%% Rafinare prin post procesare(pem) pt iv
%iau cea mai rea identificare(arx in cazul meu) si fac cu metoda err de predictie(pem)
Mvi_pem=pem(Mvi,d_id)
%validare statistica
resid(d_id,Mvi_pem)
figure
%gradul de suprapunere
compare(d_id,Mvi_pem)
%% Model cu output error (e bun!)
Moe = oe(d_id,[2,2,1])             %nb,na,nk inversate primele 2 fata de ce aveam mai sus!!
Hz = tf(Moe.B, Moe.F,dt) %fctia de transfer in discret
Hs = d2c(Hz, 'zoh') % fctia de transfer in continuu %apare zero de la 'zoh'
Moe_ss=idss(Moe)
figure;resid(d_id,Moe_ss,5)
figure;compare(d_id,Moe_ss)        

%% model cu armax da bine autocorelatia si intercorelatia(dar nu ma intereseaza intercorelatia)(E bun!)
Marmax = armax(d_id,[2,2,2,1])
Hz = tf(Marmax.B, Marmax.A,dt)  %.B si .A - proces
Hs = d2c(Hz, 'zoh')
figure;resid(d_id,Marmax)
figure;compare(d_id,Marmax)                       %Marmax = 97.23% e mai bun la resid?

%% metoda minimizarii erorii de predictie(model obtinut cu pem), intercorelatia trece, autocrelatia nu trece
Mpem = pem(d_id, 2)                                     %ordinul 2 pot spune ca trece, )))
figure
resid(d_id, Mpem)
figure
compare(d_id,Mpem)

%% metoda de minimizare a err de pred obtinut cu n4sid(trebe sa dea autocorelatia)
Mn4sid = n4sid(d_id, 1:10) % aflu ordinul sistemului
figure
resid(d_id, Mn4sid)
figure
compare(d_id,Mn4sid) %comparatie intre toate modelele


