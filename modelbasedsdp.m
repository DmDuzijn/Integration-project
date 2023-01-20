%assinging values to system parameters
Bs=1;
M_1=1;
M_2=1;
D_1=1;
D_2=1;
T=0.1 ;

%defining system matrices 
A =[1 -1*T 1*T; -(Bs/M_1)*T 1-(D_1/M_1)*T 0; (Bs/M_2)*T 0 1-(D_2/M_2)*T];
B =[0 0;(-1/M_1)*T 0;0 (-1/M_2)*T];
C= [0 1 0; 0 0 1 ];
Qy = [1 0;0 1]; %cost matrix frequency deviations
Qu = [1 0;0 1]; %cost matrix control


% starting cvx
cvx_begin sdp
variable W(3,3) symmetric %definining symmetric P matrix  
variables Y(2,3) X(2,2) %defininging dimensions of Y and X
minimize trace(transpose(C)*Qy*C*W)+trace(X) %objective
subject to
[W-eye(3) A*W+B*Y; transpose(A*W+B*Y) W] >=  0.001*eye(6); %LMI 1
W>= eye(3); %LMI 2 
[X sqrt(Qu)*Y; transpose(Y)*sqrt(Qu) W]>= 0.001*eye(5); %LMI 3
cvx_end;
%find gain matrix
K = Y*inv(W);

x0 = [0.8;0.3;0.4];legend 
%obtaining output sequences to plot open-loop system
y_o = zeros(2,41);
x_o = A*x0;
for i=2:41
    x_o = A*x_o;
    y_o(:,i)=C*x_o;
end
y_o(:,1)= C*x0;
%plotting response to random intitial condition for open loop system
t = [0:0.1:4];
plot(t,y_o,'-o')
ylabel('frequency deviation')
xlabel('time')
title('Response open loop system to random initial condition')
legend('y_1','y_2')
legend('FontSize',12)
grid

%obtaining output sequences to plot closed-loop system
y = zeros(2,41);
x = (A+B*K)*x0;
for i=2:41
    x = (A+B*K)*x;
    y(:,i)=C*x;
end
y(:,1)= C*x0;
%plotting response to random intitial condition for closed loop system
t = [0:0.1:4];
plot(t,y,'-o')
ylabel('frequency deviation')
xlabel('time')
title('Response closed loop system to random initial condition')
legend('y_1','y_2')
legend('FontSize',12)
grid





