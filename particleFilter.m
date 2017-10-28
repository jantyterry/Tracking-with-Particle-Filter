clear;

%load data
H = importdata('datafile.txt');
ydata = H(:,3);
trueposition = H(:,1);

%initialization
Qa = 0.0625;
xm1 = -10;
xm2 = 10;
Qn = 0.003906;
T = 1;

M = 1000;
E_xt = 0;
x_P = zeros(M,1);
x_P_update = zeros(M,1);
v_P = zeros(M,1);
v_P_update = zeros(M,1);
Y_t_m = zeros(M,1);
P_ytm = zeros(M,1);
w_P = (1/M)*ones(M,1);
w_P_update = (1/M)*ones(M,1);

for t = 1:numel(ydata)

%State Transition Equations
x_P_update = x_P + (v_P * T);
for i=1:M
if(x_P(i) < -20)
v_P_update(i) = 2;
elseif (x_P(i) >= -20 && x_P(i) < 0)
v_P_update(i) = v_P(i) + abs(normrnd(0,(Qa)));
elseif (x_P(i) >= 0 && x_P(i) <= 20)
v_P_update(i) = v_P(i) - abs(normrnd(0,(Qa)));
elseif (x_P(i)>20)
v_P_update(i) = -2;
end
end
x_P = x_P_update;
v_P = v_P_update;


p_1 = (1 /((sqrt(2*pi))*4));
p_2 = p_1 * exp ( (-1*((x_P - xm1).^2)) / (2 * 16) );
p_3 = p_1 * exp ( (-1*((x_P - xm2).^2)) / (2 * 16) );

%with these new updated particle locations, update the observations
%for each of these particles.
Y_t_m = p_2 + p_3;

q_1 = (1/((sqrt(2*pi))*0.003906));
P_ytm  = q_1*(exp((-(Y_t_m - ydata(t)).^2)/(2 * (2^-8)* (2^-8))));


%Using the new measurement vector yt
%the weight for each particle is updated:
w_P_update = w_P.*P_ytm;


% Normalize to form a probability distribution (i.e. sum to 1).
 nw_P_update = w_P_update ./ sum(w_P_update);

% Show distrubution of Particle with Weight before resampling
%  figure(1)
%     clf
%     plot(x_P,nw_P_update,'.k')
%     axis([-8,8,0,8.5*10^-3])
%     xlabel('Estimated particle Position')
%     ylabel('weight magnitude')
%     legend('Position of particle with weights','Location','north')
%  
 
 %Calculating Expected Values
 
 E_xt = sum(x_P .* nw_P_update);
 Ees_x(t) = E_xt;
 
 
% % %Check if sampling is necessary, and if so, resample
CV = (1/M) * sum(((M.*nw_P_update) - 1).^2);
    ESS = (M/(1 + CV));    
if (ESS < (0.5*M))  
 Q = cumsum(nw_P_update);
 T_1 = rand(M+1,1);
 T_2 = sort(T_1);
 T_2(M+1) = 1;
 a = 1; b = 1;
 while(a <= M)
     if(T_2(a)<Q(b))
         Index(a) = b;
         a = a + 1;
     else
         b = b + 1;
     end
 end
 

for a = 1:M
    
    x_P(a) = x_P(Index(a));
    v_P(a) = v_P(Index(a));
    nw_P_update(a) = 1/M;
end
end

% Show distrubution of Particle with Weight after resampling
%  figure(2)
%     clf
%     plot(x_P,nw_P_update,'.k')
%     axis([-8,8,0,8.5*10^-3])
%     xlabel('Estimated particle Position')
%     ylabel('weight magnitude')
%     legend('Position of particle with weights','Location','north')
%  

end 


% %Displaying output
figure(3)
clf
tt = 1:1109;
plot(tt,trueposition,'b')
xlabel('Time Step')
ylabel('Position')
hold on
plot(tt,Ees_x,'--r')
legend({'Trueposition','Particle Filter'},'Location','north')
hold off


figure(4)
clf
tt = 1:1109;
plot(tt,ydata,'-r')
xlabel('Time Step')
ylabel('Measurement Position')
legend('Measurement Data','Location','north')


   
