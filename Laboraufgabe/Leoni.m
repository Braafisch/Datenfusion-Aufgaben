clearvars;
in = readtable('Leoni.dat','Delimiter','semi');
y = in.Schlusskurs;
u = zeros(length(y),1);
t = in.Eintrag;
d = datenum(in.Datum,'dd.mm.yyyy');
Ts = t(2)-t(1);
R = 0.8982; Q = 1;
Ad=[1 Ts; 0 1]; Bd=[0 0]'; C=[1 0]; D=0; Gd=[Ts 1]';
%%% Init %%%
x_dach = [y(1); 0];
P_dach = 50*[1 0; 0 1];
%%% Kalman %%%
for k=1:length(y)
d_y(k) = y(k) - ( C*x_dach + D*u(k));
K = P_dach*C'*pinv(C*P_dach*C' + R);
x_tilde = x_dach + K*d_y(:,k);
P_tilde = (eye(length(Bd)) - K*C)*P_dach;
x_dach = Ad*x_tilde + Bd*u(k);
P_dach = Ad*P_tilde*Ad' + Gd*Q*Gd';
s(k) = x_tilde(1); v(k) = x_tilde(2);
K1(k) = K(1); K2(k) = K(2);
P_tilde1(k)=P_tilde(1); P_tilde2(k)=P_tilde(2);
P_tilde3(k)=P_tilde(3); P_tilde4(k)=P_tilde(4);
end
figure(1); clf;
subplot(2,1,1);
plot(t,y,'k*',t,s,'b--',...
t,s+1.5*sqrt(P_tilde1),'b-',...
t,s-1.5*sqrt(P_tilde1),'b-');
grid on; xlabel('Zeit / Tag'); ylabel('Kurs / Euro');
subplot(2,1,2);
plot(t,v,'r-*');
grid on; xlabel('Zeit / Tag');
ylabel('Kursänderungen / Euro/Tag');