clc
clear
close all

 
  % Since there is no motion in magnets so d^2y1=d^2y2=0
syms y1 y2 u1 u2 a b c d m g yc y1_0 y2_0 u1_0 u2_0   
 a = 1.65;
  b = 6.2;
  c = 2.69;
  d = 4.2;
N = 4;
m = 0.12;
 g = 9.8;
 yc = .12;
 y1_0=.02;
 y2_0=-0.02;
%  Calculation of u10 and u20
 
 u1_0= a*(b + y1_0)^4*(g*m + c/(d - y1_0 + y2_0 + yc)^4)
 u2_0=a*(g*m - c/(d - y1_0 + y2_0 + yc)^4)*(b - y2_0)^4


E= (4*c) /(((d - y1_0 + y2_0 + yc)^5))
F= (4*u1_0)/(a*(b + y1_0)^5)
G= 1       /(a*(b + y1_0)^4)
H= (4*u2_0)/(a*(b - y2_0)^5)
I= 1       /(a*(b - y2_0)^4)

% ***************************************************88
fprintf( '\n ------------------------------------------------------------------------------------------------------------\n')

 fprintf('x=Ax+Bu\n') 

 fprintf('Y=Cx+Du\n')

 
%  [0 1 0 0, -(F+E)/m  0  E/m 0 , 0 0 0 1, E/m 0 (H-E)/m 0 ]
 fprintf('A=\n')
 fprintf('   0        1     0        0 \n')
 fprintf('-(F+E)/m    0    E/m       0 \n')
 fprintf('   0        0     0        1 \n')
 fprintf('   E/m      0   (H-E)/m    0 \n')
 fprintf('B=\n')
 fprintf('  0       0   \n')
 fprintf(' G/m      0    \n')
 fprintf('  0       0    \n')
 fprintf('  0      I/m   \n')
 fprintf('\n ------------------------------------------------------------------------------------------------------------\n')
A=[0 1 0 0; -(F+E)/m  0  E/m 0 ; 0 0 0 1; E/m 0 (H-E)/m 0 ]
 B=[0 0 ; G/m 0; 0 0 ; 0 I/m]
 C=[1 0 0 0; 0 0 1 0 ]
 D=[ 0 0 ; 0 0 ]

 fprintf('\n ------------------------------------------------------------------------------------------------------------\n')
 fprintf('\n  State Space Model \n')
sys=ss(A,B,C,D)


%3 Transfer function of the open-loop system.
 fprintf('\n ------------------------------------------------------------------------------------------------------------\n')
 fprintf('\n Transfer function of the open-loop system.\n')
fprintf('           Y ( S )                   − 1       \n')
fprintf('TF ( S ) =---------- =  C ( SI − A )    B + D \n')
fprintf('           U ( S )\n')        

TF=tf(sys)
TF11=TF(1,1)
TF22=TF(2,2)
 fprintf('\n ------------------------------------------------------------------------------------------------------------\n')
 
[num1,den1] = tfdata ( TF11,'v')
[A1,B1,C1,D1] = tf2ss ( num1, den1 )
 fprintf('************************* \n')
sys1 = ss ( A1, B1, C1, D1 ) 
[num2,den2] = tfdata ( TF22,'v')
[A2,B2,C2,D2] = tf2ss ( num2, den2 )
sys2 = ss ( A2, B2, C2, D2 ) ;
fprintf('\n ------------------------------------------------------------------------------------------------------------\n')
%%4 Controllable, Observable, and Jordan canonical forms.

fprintf('---------------------------------------------------------------------------------------\n')
fprintf('3----Obtain the controllable, observable, and Jordan canonical forms.\n')
fprintf('---------------------------------------------------------------------------------------\n')
Cx=ctrb(A,B);
Rank_Cx=rank(Cx);

Cx=ctrb(A,B);
Rank_Cx=rank(Cx);
B1=B(:,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (Rank_Cx==size(A,1))
fprintf("This system is Controllable")
  CP_A =round(poly(A),6)
 n=size(A,1)
  IA1 = eye(n-1,n-1);
  Ac1=zeros(n,n);
  Bc1=zeros(n,1);
  Ac1(1:n-1,2:n)=IA1;
  Ac1(n,:)=-CP_A(:,n+1:-1:2)
  Bc1(n,1)=1
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Cx1=ctrb(A,B1)
    
    Cxc1=ctrb(Ac1,Bc1)
 
    Tc1=Cx1*inv(Cxc1)
    Cc1=C*Tc1
     Dc1=D
    if Bc1==round(inv(Tc1)*B1,4)
    
    disp('********************  Successfull ')
else 
    disp('you did wrong')
    inv(Tc1)*B1
      pause
end
% system_ctrb=ss(Ac,Bc,Cc,Dc)
% tf(system_ctrb)
else
 fprintf("This system is not Controllable")
 end
Ccanon(:,:, 1)=Cc1;
Bcanon(:,1)=Bc1
B1=B(:,2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (Rank_Cx==size(A,1))
fprintf("This system is Controllable")
  CP_A =round(poly(A),6)
 n=size(A,1)
  IA1 = eye(n-1,n-1);
  Ac1=zeros(n,n);
  Bc1=zeros(n,1);
  Ac1(1:n-1,2:n)=IA1;
  Ac1(n,:)=-CP_A (:,n+1:-1:2)
  Bc1(n,1)=1
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Cx1=ctrb(A,B1)
    
    Cxc1=ctrb(Ac1,Bc1)
 
    Tc1=Cx1*inv(Cxc1)
    Cc1=C*Tc1
     Dc1=D
    if Bc1==round(inv(Tc1)*B1,4)
    
    disp('********************  Successfull ')
else 
    disp('you did wrong')
    inv(Tc1)*B1
      pause
end
% system_ctrb=ss(Ac,Bc,Cc,Dc)
% tf(system_ctrb)
else
 fprintf("This system is not Controllable")
end
 Ccanon(:,:, 2)=Cc1
Bcanon(:,2)=Bc1


Ox=obsv(A,C)
Rank_Ox=rank(Ox);
if (Rank_Ox==size(A,1))
fprintf("This system is Observable\n")
  fprintf('Observable Canonical Form')
  Ao= Ac1'
  Bo1= Ccanon(:,:,1)'
  Bo2=Ccanon(:,:,2)'
 Co=Bcanon'
  Do=Dc1'

else
 fprintf("***** This system is not Observable")
 end

% %%%%%%%%%%%%%%
fprintf('Jordan Cnonical Form')
 [ M, R ] = eig(A)
 Aj = jordan(A)
 Bj = inv(M)*B
 Cj = C*M
 Dj = D


fprintf('---------------------------------------------------------------------------------------\n')
fprintf('4----Obtain the impulse response and the step response in Step (1) with arbitrary initial conditions.\n')
fprintf('---------------------------------------------------------------------------------------\n')
 figure
 impulse(sys1);grid
 figure
  step(sys1);grid
 figure
 impulse(sys2);grid
 figure
  step(sys2);grid


fprintf('---------------------------------------------------------------------------------------\n')
fprintf('5----Plot the Bode plot of the uncompensated system as well as the root-locus of the open-loop system\n')
fprintf('---------------------------------------------------------------------------------------\n')
 figure
 bode(TF11,'r:');grid
 figure
 rlocus(TF11,'r:');grid
 figure
 bode(TF22,'r:');grid
 figure
 rlocus(TF22,'r:');grid
% %PID Design 
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('6----Designing PID Controller\n')
fprintf('---------------------------------------------------------------------------------------\n')
  RiseTime=.25
  SettlingTime=0.005
  Overshoot=5
i=10;
while 1 
opts1 = pidtuneOptions('CrossoverFrequency',i,'PhaseMargin',60);
[ConTF11, info1] = pidtune(TF11, 'pid',opts1);
OutCsys1=feedback(ConTF11*TF11,1);
stepinfo1=stepinfo(OutCsys1);
i=i+10;
if ( stepinfo1.SettlingTime<SettlingTime)
 fprintf('Successful \n')
[ConTF11, info1] = pidtune(TF11, 'pid',opts1)
%OutCsys1=feedback(ConTF11*TF11,1);
stepinfo1
break;
end
 end
% 
%%******PID Design TF22
  i=10;

while 1 
opts2 = pidtuneOptions('CrossoverFrequency',i,'PhaseMargin',60);
[ConTF22, info2] = pidtune(TF22, 'pid',opts2);
% 
OutCsys2=feedback(ConTF22*TF22,1);
stepinfo2=stepinfo(OutCsys2);
i=i+10;
if  (stepinfo2.SettlingTime<SettlingTime)
 fprintf('Successful\n')
[ConTF22, info2] = pidtune(TF22, 'pid',opts2) 
%OutCsys2=feedback(ConTF22*TF22,1);
stepinfo2
break;
end
end

fprintf('---------------------------------------------------------------------------------------\n')
fprintf('7----Obtain the step response, square wave and sinusoidal responses.\n')
fprintf('---------------------------------------------------------------------------------------\n')
  
      tfOutCsys1=tf([ConTF11.Kd ConTF11.Kp ConTF11.Ki],[1 0])
  closelooptf1=tf(OutCsys1)
  figure
  step(OutCsys1);grid
title("Step Response of System 1 with PID  ")


   tfOutCsys2= tf([ConTF22.Kd ConTF22.Kp ConTF22.Ki],[1 0])
 closelooptf2=tf(OutCsys2)
  figure
  step(OutCsys2);grid
title("Step Response of System 2 with PID  ")


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%7 Obtain the step response, square wave and sinusoidal responses
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 7-1 sin wave TF11
   TAUg=2.5;   %Time duration 
  TFg=10;      %Spacing
  TSg=0.0001; %Time samples

% %  Generating Sinusoidal wave 
   [ Usin, Tsin ] =gensig('sin',TAUg,TFg,TSg);
%  %  Generating Squre wave
   [ Usq, Tsq   ] =gensig('squre',TAUg,TFg,TSg);
%  
%  %   Sinusoidal Wave Response of system 1 , 2 with controller
[ Ysin1, Xtsin1 ] =lsim(OutCsys1,Usin, Tsin);
[ Ysin2, Xtsin2 ] =lsim(OutCsys2,Usin, Tsin); 
%  %   Squre Wave Response of system 1 , 2  
[ Ysq1, Xtsq1   ] =lsim(OutCsys1, Usq, Tsq);
[ Ysq2, Xtsq2   ] =lsim(OutCsys2, Usq, Tsq);  
% 
% % % 
   figure
  step(OutCsys1);grid
 figure
plot( Xtsin1,Ysin1,'r','LineWidth' , 2);grid
hold on
plot( Tsin,Usin,'b','LineWidth' , 1);grid
hold off   
grid
title("Sinusoidal Wave Response of system 1 with PID")

% 
% % 7-2 squre TF11
% 
% % %
figure
 plot( Xtsq1, Ysq1,'r');grid
 hold on
 plot( Tsq, Usq,'b-');grid
  hold off   
 grid
 
title("Squre Wave Response of system 1 with PID")
% % %%%%%%%%%%%%%%%%%%sin wave TF22
% 
%  
% % % 
   figure
  step(OutCsys2);grid
figure
plot( Xtsin2,Ysin2,'r','LineWidth' , 2);grid
hold on
plot( Tsin,Usin,'b-','LineWidth' ,  1);grid
  hold off  
 grid
title("Sinusoidal Wave Response of system   2 with PID")
% 
% % %%%%%%%%%%%%%%%%%%squre TF22
% 
%  
% 
% % % 
figure
 plot( Xtsq2, Ysq2,'r');grid   
 hold on
 plot( Tsq, Usq,'b-');grid
 hold off  
 grid
title("Squre Wave Response of system   2 with PID") 

fprintf('---------------------------------------------------------------------------------------\n')
fprintf('8----Plot the control input signal for each set point input .\n')
fprintf('---------------------------------------------------------------------------------------\n')
   figure
   plot(Tsin,Usin-Ysin1) ;grid
title("Control Input Signal for Set Point Sinusoidal Wave input1")
    figure
   plot(Tsin,Usin-Ysin2) ;grid
 title("Control Input Signal for Set Point Sinusoidal Wave input2")
   figure
    plot(Tsq,Usq-Ysq1) ;grid
title("Control Input Signal for Set Point squre Wave input1")
     figure
    plot(Tsq,Usq-Ysq2);grid
title("Control Input Signal for Set Point squre Wave input2")
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% 9 %% Examine the robustness of the design by introducing noise and parameter variations or uncertainty in the system.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('9----Examine the robustness of the design by introducing noise and parameter variations or uncertainty in the system.\n')
fprintf('---------------------------------------------------------------------------------------\n')
  figure
  rlocus(OutCsys1,'r-');grid
    figure
  rlocus(OutCsys2,'r-');grid
  
    % Noise and Disturbance
 % respone step Noise and Disturbance Controller 1
contr1=tf([ConTF11.Kd ConTF11.Kp ConTF11.Ki],[1 0]);
loop_contr1 = TF11 * contr1;
loop_distr1 = 1 + loop_contr1;
resp1 = TF11 / loop_distr1;
figure
step(resp1);grid
title("Step respone to Noise and Disturbance System 1")



 % respone step Noise and Disturbance Controller 2
contr2=tf([ConTF22.Kd ConTF22.Kp ConTF22.Ki],[1 0]);
loop_contr2 = TF22 * contr2;
loop_distr2 = 1 + loop_contr2;
resp2 = TF22 / loop_distr2;
figure
step(resp2);grid
title("Step respone to Noise and Disturbance System 2")



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   10 %% Design a full state feedback control to meet the design   %
%         specification indicated in Step(6).                       %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('10----Design a full state feedback control to meet the design specification indicated.\n')
fprintf('---------------------------------------------------------------------------------------\n')
    SettlingTime1=.005
    SettlingTime2=0.05
    z1=1
    z2=0.707
   Wn1=4/(SettlingTime1*z1)
   Wn2=4/(SettlingTime2*z2)
  Pf=zeros(1,4);
  Pf(1,1:2)= roots([1 2*z1*Wn1 Wn1^2]).';%s2+2zwns+wn2
  Pf(1,3:4)= roots([1 2*z2*Wn2 Wn2^2]).'

   K=place(A,B ,Pf) %Closed-loop pole assignment using state feedback

  Af= A- B*K
pause
%x'=Afx+Bu  y=-kx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%11 Plot the control input signal for each set point input in Step (7).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('11----Plot the control input signal for each set point input.\n')
fprintf('---------------------------------------------------------------------------------------\n')
 sysnew=ss(Af ,B ,C,D)
 Kdc_c=dcgain(sysnew)
 Kr_c=inv(Kdc_c)

%   Scaling the system 
  sysscaled=ss(Af,B*Kr_c,C,D);
%   Making new transfer function 
   Tfnew=tf(sysscaled);
%  Step response to the system   
   figure
   step(sysscaled);grid
title("Step response of State Feedback system");
%  input 1 and output 1
sys1new=ss(sysscaled.A,sysscaled.B(:,1),sysscaled.C(1,:),sysscaled.D(1,1))
%  input 2 and output 2
sys2new=ss(sysscaled.A,sysscaled.B(:,2),sysscaled.C(2,:),sysscaled.D(2,2))

%   Sinusoidal Wave Response of system 1
[ Ynewsin1, Xnewtsin1 ] =lsim(sys1new,Usin, Tsin);

%   Sinusoidal Wave Response of system 2
[ Ynewsin2, Xnewtsin2 ] =lsim(sys2new,Usin, Tsin);
%   Squre Wave Response of system 1
[ Ynewsq1, Xnewtsq1 ] =lsim(sys1new,Usq, Tsq);
%   Squre Wave Response of system 2
[ Ynewsq2, Xnewtsq2 ] =lsim(sys2new,Usq, Tsq);
%  Drawnin and Comparing  Sinusoidal Wave between input and output system 1 
figure
  plot( Tsin,Ynewsin1,'r','LineWidth' , 2);grid
  hold on   
  plot( Tsin,Usin,'b','LineWidth' , 1);grid
  hold off  
 grid
title("Sinusoidal Wave of State Feedback system 1 ");


%  Drawnin and Comparing  Sinusoidal Wave between input and output system 2
figure
   plot( Tsin,Ynewsin2,'r','LineWidth' , 2);grid
   hold on
   plot( Tsin,Usin,'b','LineWidth' , 1);grid
hold off  
 grid
title(" Sinusoidal Wave of State Feedback system 2");

%  Drawnin and Comparing  Squre Wave between input and output system 1
figure
plot( Xnewtsq1,Ynewsq1,'r','LineWidth' , 2);grid
 hold on
plot( Tsq,Usq,'b-','LineWidth' ,  1);grid
hold off  
 grid
title("Squre Wave of State Feedback system 1 ");

  %  Drawnin and Comparing  Squre Wave between input and output system 2

figure
plot( Xnewtsq2,Ynewsq2,'r','LineWidth' , 2);grid
hold on
plot( Tsq,Usq,'b-','LineWidth' ,  1);grid
  hold off  
 grid
title("Squre Wave of State Feedback system 2 ");


 figure
 rlocus(sys1new,'r:');grid
 figure
 rlocus(sys2new,'r:');grid
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %12Plot the control input signal for each set point input in Step (11)
 %and compare the result with those in Step (8).
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('11----Plot the control input signal for each set point input in Step (11)and compare the result with those in Step (8).\n')
fprintf('---------------------------------------------------------------------------------------\n')
    figure
   plot(Tsin,Usin-Ynewsin1) ;grid
   title("Sinusoidal Input Control response system 1")
  figure
   plot(Tsin,Usin-Ynewsin2);grid 
  title("Sinusoidal Input Control response system 2")
  
   figure
    plot(Tsq,Usq-Ynewsq1) ;grid
    
  title("Squre Input Control response system 1")
  figure
    plot(Tsq,Usq-Ynewsq2);grid
    title("Squre Input Control response system 2")
    
%  
%  
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  13 Design a full-order and a reduced-order observer 
%  %and obtain step and sinusoidal responses.
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('13----Design a full-order and a reduced-order observer and obtain step and sinusoidal responses.\n')
fprintf('---------------------------------------------------------------------------------------\n')
  fprintf('13-1---Design a full-order obtain step and sinusoidal responses.') 
% %%13-1 Full_order observer
% % %%%%%%%%
    SettlingTime1=SettlingTime1;
    SettlingTime2=SettlingTime2;
  f_time=4; 
  z1_obs=1;
  z2_obs=0.707;
% 
% 
SettlingTime_obs1=SettlingTime1/f_time
SettlingTime_obs2=SettlingTime2/f_time
   Wn1_obs=4/(SettlingTime_obs1*z1_obs)
   Wn2_obs=4/(SettlingTime_obs2*z2_obs)
  P_obs=zeros(1,4);
  P_obs(1,1:2)= roots([1 2*z1_obs*Wn1_obs Wn1_obs^2])'
  P_obs(1,3:4)= roots([1 2*z2_obs*Wn2_obs Wn2_obs^2])'
  G =place(A',C' ,P_obs)'
% 
% 
    A_est=A-(G*C)-(B*K)
    B_est=G
    C_est=-K
    D_est=0
   sys_obs=ss(A_est,B_est,C_est,D_est)
   Tf_sys_obs=tf(sys_obs)
   Kdc_o=dcgain(sys_obs);
   Kr_o=inv(Kdc_o)
     sys_obs_scaled=ss(A_est,B_est*Kr_o,C_est,D_est)
    figure
   step( sys_obs_scaled );grid
%  
  sys1_obs_con=ss(sys_obs_scaled.A,sys_obs_scaled.B(:,1),sys_obs_scaled.C(1,:),sys_obs_scaled.D(1,1))
% 
  sys2_obs_con=ss(sys_obs_scaled.A,sys_obs_scaled.B(:,2),sys_obs_scaled.C(2,:),sys_obs_scaled.D(2,2))
% 
% Sinusoidal reponse to full oreder observer system 1 
  [ Y_ocsin1, X_octsin1 ] =lsim(sys1_obs_con,Usin , Tsin);
  [ Y_ossq1, Xostsq1    ] =lsim(sys1_obs_con,Usq, Tsq);
% Sinusoidal reponse to full oreder observer system 2 
  [ Y_ocsin2, X_octsin2 ] =lsim(sys2_obs_con,Usin , Tsin);
  [ Y_ossq2, Xostsq2    ] =lsim(sys2_obs_con,Usq, Tsq);

figure
  plot( Tsin,Y_ocsin1,'r','LineWidth' , 2);grid
  hold on   
  plot( Tsin,Usin,'b','LineWidth' , 1);grid
  hold off  
 grid
 title(" Sinusoidal Wave reponse of Full Order Observer system 1");

 figure
  plot( Tsin,Y_ocsin2,'r','LineWidth' , 2);grid
  hold on   
  plot( Tsin,Usin,'b','LineWidth' , 1);grid
  hold off  
 grid
title(" Sinusoidal Wave reponse of Full Order Observer system 2");


figure
plot( Xostsq1,Y_ossq1,'r','LineWidth' , 2) 
 hold on
plot( Tsq,Usq,'b-','LineWidth' ,  1) 
hold off
grid
title(" Squre Wave reponse of Full Order Observer system 1");
  
figure
plot( Xostsq2,Y_ossq2,'r','LineWidth' , 2);grid
hold on
plot( Tsq,Usq,'b-','LineWidth' ,  1);grid
  hold off
grid
 title(" Squre Wave reponse of Full Order Observer system 2");
fprintf('*******************************************************************.\n')
fprintf('13-2---Design a reduced-order obtain step and sinusoidal responses.\n') 
% 
% % %13-2 Reduced Observer   
% %  
% % %%alpha_e(s)=
    SettlingTime1_robsv=SettlingTime1/f_time
    SettlingTime2_robsv=SettlingTime2/f_time
    z1=1
    z2=0.707
   Wn1_r=4/(SettlingTime1_robsv*z1)
   Wn2_r=4/(SettlingTime2_robsv*z2)
  Pfnew=zeros(1,2);
n=size(A,1)
p=size(B,2)
q=size(C,1)
P=[1 0 0 0; 0 0 1 0;0 1 0 0 ;0 0 0 1]
pinverse=inv(P)
A_new=P*A*inv(P)
B_new=P*B
C_new=C*inv(P)
K_new=K*inv(P)

A11=A_new(1:q,1:q)
A12=A_new(1:q,q+1:n)
A21=A_new(q+1:n,1:q)
A22=A_new(q+1:n,q+1:n)
B1=B_new(1:q,:)
B2=B_new(q+1:n,:)
K1=K_new(:,1:q)
K2=K_new(:,q+1:n)
  Pfnew(1,1:2) =roots([1 2*z1*Wn1_r Wn1_r^2])%s2+2zwns+wn2
%    Pfnew(1,3:4)=roots([1 2*z2*Wn2_r Wn2_r^2])'
% 
 Gnew=place(A22',A12' ,Pfnew )' %Closed-loop pole assignment using state feedback
A_cre=A22-(Gnew*A12)-(B2*K2)+(K2*B1*K2)
B_cre=(A22*Gnew)-(Gnew*A12*Gnew)-(Gnew*A11)+A21-(B2*(K1+(K2*Gnew)))+((Gnew*B1)*(K1+(K2*Gnew)))
C_cre=-K2
D_cre=-(K1+(K2*Gnew))
  sys_sys_cre = ss ( A_cre, B_cre, C_cre, D_cre ) ;
   Tf_sys_cre=tf(sys_sys_cre)
   Kdc_cre=dcgain(Tf_sys_cre);
%     Kr_cre=inv(Kdc_cre)
      Kr_cre=1
     
     sys_cre_scaled=ss(A_cre,B_cre*Kr_cre,C_cre,D_cre)
    figure
   step( sys_cre_scaled ) ;grid
% % 
  sys1_cre_con=ss(sys_cre_scaled.A,sys_cre_scaled.B(:,1),sys_cre_scaled.C(1,:),sys_cre_scaled.D(1,1))
% 
  sys2_cre_con=ss(sys_cre_scaled.A,sys_cre_scaled.B(:,2),sys_cre_scaled.C(2,:),sys_cre_scaled.D(2,2))
% 
% Sinusoidal reponse to full oreder observer system 1 
  [ Y_cresin1, X_cretsin1 ] =lsim(sys1_cre_con,Usin , Tsin);
  [ Y_cresq1, Xcretsq1    ] =lsim(sys1_cre_con,Usq, Tsq);
% Sinusoidal reponse to full oreder observer system 2 
  [ Y_cresin2, X_cretsin2 ] =lsim(sys2_cre_con,Usin , Tsin);
  [ Y_cresq2, Xcretsq2    ] =lsim(sys2_cre_con,Usq, Tsq);

figure
  plot( Tsin,Y_cresin1,'r','LineWidth' , 2);grid
  hold on   
  plot( Tsin,Usin,'b','LineWidth' , 1);grid
  hold off  
 grid
 title(" Sinusoidal Wave reponse of reduced Order Observer system 1");

 figure
  plot( Tsin,Y_cresin2,'r','LineWidth' , 2);grid
  hold on   
  plot( Tsin,Usin,'b','LineWidth' , 1);grid
  hold off  
 grid
title(" Sinusoidal Wave reponse of reduced Order Observer system 2");


figure
plot( Xcretsq1,Y_cresq1,'r','LineWidth' , 2);grid
 hold on
plot( Tsq,Usq,'b-','LineWidth' ,  1);grid
hold off  
 grid
 title(" Squre Wave reponse of reduced Order Observer system 1");
  
figure
plot( Xcretsq2,Y_cresq2,'r','LineWidth' , 2);grid
hold on
plot( Tsq,Usq,'b-','LineWidth' ,  1);grid
  hold off  
 grid
 title(" Squre Wave reponse of reduced Order Observer system 2");% % 
% % % 
% % % %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %  %%%%14 Find the transfer function of the observer and controller.
% % % %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 fprintf('---------------------------------------------------------------------------------------\n')
fprintf('14----Find the transfer function of the observer and controller.\n')
fprintf('---------------------------------------------------------------------------------------\n')
sympref ( 'FloatingPointOutput' , true ) ;
syms s   
   


 TF_cfo =tf(sys_obs_scaled)
TF_cre =tf(sys_cre_scaled)
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('----Saving diagrams .\n')
fprintf('---------------------------------------------------------------------------------------\n')
  for i=1:45 
  figurename="figure"+int2str(i);
  saveas(figure(i), figurename, 'jpg');
 
  end
close all
fprintf('------------------Finish---------------- .\n')
