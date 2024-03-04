t=Pirvan(:,1)
u=Pirvan(:,2)  %intrarea
y2=Pirvan(:,4) %a doua iesire
plot(t,u,t,y2), title('Pentru a doua iesire')

%%
%de pe semnalul de iesire
i1=412;
i2=432;
%de pe semnalul de intrare
i3=405;
i4=427;

%factorul de proporționalitate în regim staționar
k = mean(y2)/mean(u)
%intarziere intre semnalul de intrare si semnalul de iesire(defazajul?)
dt = t(i1)-t(i3)
%amplificarea la rezonanta
Mr = (y2(i1)-y2(i2))/(u(i3)-u(i4))/k
%perioada la rezonanta pentru semnalul de intrare
Tr = 2*(t(i4)-t(i3))
%pulsatia la rezonanta
wr = 2*pi/Tr
%factor de amortizare
zeta =sqrt((Mr-sqrt(Mr^2-1))/2/Mr)
%pulsatia oscilatiilor
wn = wr/sqrt(1-2*zeta^2)
 
%defazajul la rezonanta
phr = (t(i3)-t(i1))*wr %in radiani
phr = (t(i3)-t(i1))*wr*180/pi %in grade
%constanta de timp a zeroului
Tz = tan(phr + atan(sqrt(1 - 2*zeta^2)/zeta))/wr
 
num = k*wn^2*[Tz, 1]
den = [1, 2*zeta*wn, wn^2]
%param Markov
g1 = k*Tz*wn^2
g2 = k*wn^2-k*Tz*wn^3*2*zeta
%spatiul starilor
A = [0, 1; -wn^2, -2*zeta*wn]
B = [g1 ;g2]
C = [1 0]
D = [0]
dt = t(2)-t(1);
yc = lsim(A,B,C,D, u,t, [y2(1), (y2(2)-y2(1))/dt - g1*u(1)]);
plot(t, [y2 yc])

J1 = norm(y2-yc)/sqrt(length(y2))
eMPN1 = norm(y2-yc)/norm(y2-mean(y2))*100

%% Identificare parametrica(al doilea cel mai prost -> fac rafinare)
dt2 = t(2)-t(1) %pas de achizitie
d_id = iddata(y2,u,dt)
%model calculat cu arx                                                      (identificare prin arx)
Marx = arx(d_id,[2,2,1])
figure;resid(d_id, Marx, 5) 
figure;compare(d_id, Marx) 
Hz2 = tf(Marx.B, Marx.A,dt2) %fctia de transfer in discret
Hs2 = d2c(Hz2, 'zoh')% discrete to continous

%% Rafinarea lui arx prin metoda minimizarii erorii de predictie
Marx_pem=pem(Marx,d_id)
resid(d_id,Marx_pem,5)
figure
compare(d_id,Marx_pem)
%% Identif cu variab instrumentale (CEL MAI PROST)
Mvi = iv4(d_id, [2,2,1])
figure;resid(d_id, Mvi, 5)
figure;compare(d_id,Mvi)

%% rafinare cu pem pt iv  (eroarea s a miscorat putin)
Mvi_pem=pem(Mvi,d_id)
%validare statistica
resid(d_id,Mvi_pem,5)
figure
compare(d_id,Mvi_pem)

%% Model cu output error (BUN!) trece intercor
Moe = oe(d_id,[2,2,1])
Hz = tf(Moe.B, Moe.F,dt) %fctia de transfer in discret
Hs = d2c(Hz) 
figure;resid(d_id,Moe,5)
figure;compare(d_id,Moe)

%% model cu armax (BUN!) trece autocor
Marmax = armax(d_id,[2,2,2,1])
Hz = tf(Marmax.B, Marmax.A,dt)  
Hs = d2c(Hz, 'zoh')
figure;resid(d_id,Marmax)
figure;compare(d_id,Marmax)

%% metoda minimizarii erorii de predictie(model obtinut cu pem), intercorelatia trece, autocrelatia nu trece
%trece autocorelatia
Mpem = pem(d_id, 2)                               %ordinul 2 pot spune ca trece, )))
figure
resid(d_id, Mpem)
figure
compare(d_id,Mpem)

%% metoda de minimizare a err de pred obtinut cu n4sid(trebe sa dea autocorelatia)
Mn4sid = n4sid(d_id, 1:10)
figure
resid(d_id, Mn4sid,5)
figure
compare(d_id,Mn4sid)

