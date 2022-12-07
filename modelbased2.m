Bs=1;
M_1=1;
M_2=1;
D_1=1;
D_2=1;
T=0.1 ;

A =[1 -1*T 1*T; -(Bs/M_1)*T 1-(D_1/M_1)*T 0; (Bs/M_2)*T 0 1-(D_2/M_2)*T];
B =[0 0;(-1/M_1)*T 0;0 (-1/M_2)*T]
C= [0 1 0; 0 0 1 ];
z= [0;1;0]
Qy = eye(2);
Qu = eye(2);

% starting cvx
cvx_begin sdp
variable P(3,3) symmetric %definining symmetric P matrix  
variables Y(2,3) X(2,2) %defininging dimensions of Y and X
minimize trace(transpose(C)*Qy*C*P)+trace(X) %objective
subject to
[P-eye(3) A*P+B*Y; transpose(A*P+B*Y) P] >=  0.001*eye(6); %LMI 1
P>= eye(3); %LMI 2 
[X sqrt(Qu)*Y; transpose(Y)*sqrt(Qu) P]>= 0.001*eye(5); %LMI 3
cvx_end;

K = Y*inv(P);

A_cl=A+B*K;
B_cl=eye(3);
cl_sys=ss(A_cl,B_cl,C,0);
step(cl_sys,10)