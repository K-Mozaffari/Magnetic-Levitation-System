clear
clc
% Since there is no motion in magnets so d^2y1=d^2y2=0 
syms y1 y2 u1 u2 a b c d m g yc y1_0 y2_0 u1_0 u2_0 mg1 mg2 mg1_0 mg2_0  md2y1 md2y2 Formula_u1_0  Formula_u2_0 
N = 4;
y12 = yc + y2 - y1;
y12_0 = yc + y2_0 - y1_0;
Fu11 = u1 /(a * (y1 + b)^N);
Fu22 = u2 /(a * (-y2 + b)^N);
Fm12 = c / (y12 + d)^N;
mg1  = u1  / (a * (y1  + b)^N)-(c / (y12  + d)^N);
mg2  = u2  / (a * (-y2  + b)^N)+(c / (y12  + d)^N) ;
mg1_0 = u1_0 / (a * (y1_0 + b)^N)-(c / (y12_0 + d)^N);
mg2_0 = u2_0 / (a * (-y2_0 + b)^N)+c / (y12_0 + d)^N ;
Formula_u1_0  = - m*g   + (u1_0  / (a * (y1_0   + b)^N)-(c / (y12_0  + d)^N));
Formula_u2_0  = - m*g   + (u2_0  / (a * (-y2_0  + b)^N)+(c / (y12_0  + d)^N));
md2y1 =( Fu11 - Fm12 - mg1_0) ;
md2y2 = (Fu22 + Fm12- mg2_0)   ;
% Soluotions of linearization by using taylor series expanstion
    F1 = simplify(taylor(md2y1, [y1, y2, u1], [y1_0, y2_0, u1_0], 'Order',2));
    F2 = simplify(taylor(md2y2, [y1, y2, u2], [y1_0, y2_0, u2_0], 'Order',2));
 fprintf('my"1=\n')
pretty( F1)
 fprintf('\n ************************************************************************************************************\n')
fprintf('my"2=\n')
pretty( F2)
fprintf('\n ------------------------------------------------------------------------------------------------------------\n')
fprintf('u1_0= \n')
pretty (solve (Formula_u1_0 ,u1_0 ))
fprintf('\n ------------------------------------------------------------------------------------------------------------\n')
fprintf('u2_0= \n')
pretty (solve (Formula_u2_0 ,u2_0 ))



