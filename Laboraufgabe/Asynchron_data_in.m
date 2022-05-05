clearvars;
in = readtable('Asynchron_data_in.dat','Delimiter','space');
t = in.time;
y(:,1)=in.position;
y(:,2)=in.acceleration;
u = zeros(1,length(y));

R = 0.8982; Q = 1;
%%% Init %%%
x_tilde = [y(1,1); 0; 0];
P_tilde = 50*eye(3);
%%% Kalman %%%
for k=2:length(y)
    Ts = t(k)-t(k-1);
    Ad = [1 Ts Ts^2/2; 0 1 Ts; 0 0 1];
    Bd = [0; 0; 0];
    D = 0;
    Gd = [Ts^2/2; Ts; 1];

    if isnan(y(k,1))
        C = [0 0 1];
        y_n = y(k,2);
    else
        C = [1 0 0];
        y_n = y(k,1);
    end
    x_dach = Ad*x_tilde + Bd.*u(k-1);
    P_dach = Ad*P_tilde*Ad' + Gd*Q*Gd';

    d_y(k) = y_n - ( C*x_dach + D*u(k));
    K = P_dach*C'*pinv(C*P_dach*C' + R);

    x_tilde = x_dach + K*d_y(:,k);
    P_tilde = (eye(length(Bd)) - K*C)*P_dach;
    s(k) = x_tilde(1); v(k) = x_tilde(2); a(k) = x_tilde(3);
    Ps(k) = P_tilde(1); Pv(k) = P_tilde(5); Pa(k) = P_tilde(9);
end
h2 = figure(2); clf;
%%% Weg (Ein- und Ausgangssignal) %%%
subplot(3,1,1);
plot(t,y(:,1),'ko',t,s,'r-*',t,s+3*sqrt(Ps),'b-',t,s-3*sqrt(Ps),'b-'); grid on;
xlabel('Zeit / s'); ylabel('Position s(t) / m');
%%% Geschwindigkeit (Aussgangssignal) %%%
subplot(3,1,2);
plot(t,v,'r-*',t,v+3*sqrt(Pv),'b-',t,v-3*sqrt(Pv),'b-'); grid on;
xlabel('Zeit / s'); ylabel('Geschw. m/s');
%%% Beschleunigung (Aussgangssignal) %%%
subplot(3,1,3);
plot(t,y(:,2),'ko',t,a,'r-*',t,a+3*sqrt(Pa),'b-',t,a-3*sqrt(Pa),'b-'); grid on;
xlabel('Zeit / s'); ylabel('Beschl. m/s^2');