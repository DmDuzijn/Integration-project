%assinging the parameters to a value of 1.
Bs=1;
M_1=1;
M_2=1;
D_1=1;
D_2=1;
%choosing a time interval 0.1 seconds
T=0.1 ;

%Defining the open-loop A,B and C matrices
A =[1 -1*T 1*T; -(Bs/M_1)*T 1-(D_1/M_1)*T 0; (Bs/M_2)*T 0 1-(D_2/M_2)*T];
B =[0 0;(-1/M_1)*T 0;0 (-1/M_2)*T];
C=[0 1 0;0 0 1];
%defining the cost matrices
Qx =[0 0 0;0 1 0; 0 0 1]; 
Qy= eye(2); %only takes the outputs into account
R = eye(2); %assigns costs for control

%creating the data-sequences, note that, as stated in the report, the
%minimimal length of the data-sequences should be 11. Here, data-sequences
%of 100 are created to be able to plot the open-loop system on a larger
%interval.
u1 = rand(2,100); %a random input sequence is used
x = zeros(3,101);
y=zeros(2,101);
x0 = [0.8;0.30;0.38];
u0 = rand(2,1);
x(1:3,1)=A*x0+B*u0; %determining the first value of the state sample matrix
y(1:2,1)=C*x0; %determing the first value of the output matrix
for k=1:100 %this for loop determines the ouput value for each iteration to construct an output matrix
    x(1:3,k+1)=A*x(1:3,k)+B*u1(:,k);
    y(1:2,k+1)=C*x(1:3,k);
end

%extracting the appropiate length of the data-sequence matrices to obtain
%data-sequences which are not longer than neccessary. 

X0 = [x0 x(:,1:10)];
Y = C*X0
u = [u0 u1(:,1:10)];
X1 =[B A]*[u;X0]; %defining X1
Yn = Y(:,1:3); %extraction of first three colums to be able to express output matrix C in terms of data
Xn = X0(:,1:3); %extraction of first three colums to be able to express output matrix C in terms of data
NewC=Yn*inv(Xn) %Output matrix C could now be obtained with this formula.

d = det(Xn) %checking if Xn is invertible
if d == 0
    display('NOT INVERTIBLE')
end 

%Starting the CVX-Toolbox
cvx_begin sdp
variable Q(11,3) %Defining the dimensions of the decision variables
variable X(2,2) %Defining the dimensions of the decision variables
minimize trace(transpose(NewC)*Qy*NewC*X0*Q)+trace(X) %objective function
[X0*Q-eye(3) X1*Q;transpose(Q)*transpose(X1) X0*Q] >= 0.01*eye(6); %LMI 1
[X R^0.5*u*Q; transpose(Q)*transpose(u)*R^0.5 X0*Q ] >= 0.01*eye(5); %LMI 2

cvx_end

K = u*Q*inv(X0*Q); %retrieving gain matrix
Gk=Q*inv(X0*Q); %retrieving Gk for own understanding
KI=[u;X0]*Gk;   %varifiying equation 4.13 for own understanding


%closed loop system 1, this will be used to compare different weighting
%matrices
A_cl = X1*Q*inv(X0*Q);
B_cl = -B*K;
C_cl = Yn*inv(Xn);
D_cl = 0;
cl_sys = ss(A_cl,B_cl,C_cl,D_cl,0.1);

%closed loop system 2 and 3, this will be used to compare different weighting
%matrices
[L, I, J] = dlqr(A,B,Qx,[0.1 0;0 0.1]);
cl2_sys=ss(A-B*L,B*L,C,0,0.1);
[L2, I2, J2] = dlqr(A,B,Qx,[0.01 0;0 0.01]);
cl2_sys=ss(A-B*L,B*L,C,0,0.1);

%plotting response to random input sequence for open loop system
t= [0:0.1:4];
nexttile
plot(t,u1(:,1:41),'-o')
ylabel('Input')
xlabel('Time')
title('Input data-sequences open-loop system')
legend('u_1','u_2')
legend('FontSize',12)
legend('location','northeastoutside')
grid
nexttile
plot(t,y(:,1:41),'-o')
ylabel('Output')
xlabel('Time')
title('Output data-sequences open-loop system')
legend('y_1','y_2')
legend('FontSize',12)
legend('location','northeastoutside')
grid

%plotting response to random intitial condition for closed loop system
y_cl = zeros(2,41);
x = X1*Gk*x0;
for i=2:41
    x = X1*Gk*x;
    y_cl(:,i)=C*x;
end
y_cl(:,1)= C*x0;
t = [0:0.1:4];
plot(t,y_cl,'-o')
ylabel('Frequency deviation')
xlabel('Time')
title('Response closed loop system to random initial condition')
legend('y_1','y_2')
legend('FontSize',12)
grid

%plotting different weighting matrices
y_cl2 = zeros(2,41);
x = (A-B*L)*x0;
for i=2:41
    x = (A-B*L)*x;
    y_cl2(:,i)=C*x;
end
y_cl2(:,1)= C*x0;

y_cl3 = zeros(2,41);
x = (A-B*L2)*x0;
for i=2:41
    x = (A-B*L2)*x;
    y_cl3(:,i)=C*x;
end
y_cl3(:,1)= C*x0;

%plotting reaction of closed-loop systems desgined with different weighting
%matrices
t = [0:0.1:4];
plot(t,y_cl(1,:),'-o')
hold on
plot(t,y_cl2(1,:),'-o')
hold on
plot(t,y_cl3(1,:),'-o')
title('Response different systems to random initial condition')
ylabel('Frequency deviation bus 1')
xlabel('Time')
legend('sys1','sys2','sys3')
grid
hold off
nexttile
plot(t,y_cl(2,:),'-o')
hold on
plot(t,y_cl2(2,:),'-o')
hold on
plot(t,y_cl3(2,:),'-o')
title('Response different systems to random initial condition')
ylabel('Frequency deviation bus 2')
xlabel('Time')
legend('sys1','sys2','sys3')
grid
hold off







