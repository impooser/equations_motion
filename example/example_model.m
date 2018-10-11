clear
syms m_V m_L1 m_L2 m_L3 I_V I_L2 I_L3 L1 L2 g_0
syms theta_1(t) theta_2(t) beta(t) x(t)

%% Lagrangian

% X = x;
% Y = 0;
% 
% T1 = simplify(1/2*m_L1*(diff(X,t)^2+diff(Y,t)^2));
% U1 = simplify(m_L1*g_0*Y);
% 
% X = x-L1/2*cos(theta_1);
% Y = L1/2*sin(theta_1);
% 
% T2 = simplify(1/2*m_L2*(diff(X,t)^2+diff(Y,t)^2) + 1/2*I_L2*diff(theta_1,t)^2);
% U2 = simplify(m_L2*Y*g_0);
% 
% X = x - L1*cos(theta_1) - L2/2*cos(theta_2);
% Y = L1*sin(theta_1) + L2/2*sin(theta_2);
% 
% T3 = simplify(1/2*m_L3*(diff(X,t)^2+diff(Y,t)^2) + 1/2*I_L3*diff(theta_2,t)^2);
% U3 = simplify(m_L3*Y*g_0);
% 
% X = x - L1*cos(theta_1) - L2*cos(theta_2);
% Y = L1*sin(theta_1) + L2*sin(theta_2);
% 
% T4 = simplify(1/2*m_V*(diff(X,t)^2+diff(Y,t)^2) + 1/2*I_V*diff(beta,t)^2);
% U4 = simplify(m_V*Y*g_0);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% T = T1 + T2 + T3 + T4;
% U = U1 + U2 + U3 + U4;
% 
% L = T - U;

%% Derivatives of Lagrangian
% 
% clear
% 
% syms m_V m_L1 m_L2 m_L3 I_V I_L2 I_L3 L1 L2 g_0
% syms theta_1 theta_2 beta x dottheta_1 dottheta_2 dotbeta dotx
% 
% L = (m_L3*((dotx + L1*sin(theta_1)*dottheta_1 + (L2*sin(theta_2)*dottheta_2)/2)^2 + (L1*cos(theta_1)*dottheta_1 + (L2*cos(theta_2)*dottheta_2)/2)^2))/2 + (m_V*((dotx + L1*sin(theta_1)*dottheta_1 + L2*sin(theta_2)*dottheta_2)^2 + (L1*cos(theta_1)*dottheta_1 + L2*cos(theta_2)*dottheta_2)^2))/2 + (I_V*dotbeta^2)/2 + (I_L2*dottheta_1^2)/2 + (I_L3*dottheta_2^2)/2 + (m_L1*dotx^2)/2 + (m_L2*dotx^2)/2 + (L1^2*m_L2*dottheta_1^2)/8 - g_0*m_L3*(L1*sin(theta_1) + (L2*sin(theta_2))/2) - g_0*m_V*(L1*sin(theta_1) + L2*sin(theta_2)) - (L1*g_0*m_L2*sin(theta_1))/2 + (L1*m_L2*sin(theta_1)*dottheta_1*dotx)/2;
% 
% dL_dotx = simplify(diff(L,dotx));
% dL_x = simplify(diff(L,x));
% 
% dL_dottheta1 = simplify(diff(L,dottheta_1));
% dL_theta1 = simplify(diff(L,theta_1));
% 
% dL_dottheta2 = simplify(diff(L,dottheta_2));
% dL_theta2 = simplify(diff(L,theta_2));
% 
% dL_dotbeta = simplify(diff(L,dotbeta));
% dL_beta = simplify(diff(L,beta));

%% Equations of motion

%clear
% 
% syms theta_1(t) theta_2(t) beta(t) x(t) dottheta_1(t) dottheta_2(t) dotbeta(t) dotx(t)
% 
% ddtdL_dotx = simplify(diff(dotx*m_L1 + dotx*m_L2 + (m_L3*(2*dotx + 2*L1*dottheta_1*sin(theta_1) + L2*dottheta_2*sin(theta_2)))/2 + (m_V*(2*dotx + 2*L1*dottheta_1*sin(theta_1) + 2*L2*dottheta_2*sin(theta_2)))/2 + (L1*dottheta_1*m_L2*sin(theta_1))/2 , t));
% ddtdL_dottheta1 = simplify(diff(I_L2*dottheta_1 + (m_L3*(2*L1*sin(theta_1)*(dotx + L1*dottheta_1*sin(theta_1) + (L2*dottheta_2*sin(theta_2))/2) + 2*L1*cos(theta_1)*(L1*dottheta_1*cos(theta_1) + (L2*dottheta_2*cos(theta_2))/2)))/2 + (m_V*(2*L1*sin(theta_1)*(dotx + L1*dottheta_1*sin(theta_1) + L2*dottheta_2*sin(theta_2)) + 2*L1*cos(theta_1)*(L1*dottheta_1*cos(theta_1) + L2*dottheta_2*cos(theta_2))))/2 + (L1^2*dottheta_1*m_L2)/4 + (L1*dotx*m_L2*sin(theta_1))/2 , t));
% ddtdL_dottheta2 = simplify(diff(I_L3*dottheta_2 + (m_L3*(L2*sin(theta_2)*(dotx + L1*dottheta_1*sin(theta_1) + (L2*dottheta_2*sin(theta_2))/2) + L2*cos(theta_2)*(L1*dottheta_1*cos(theta_1) + (L2*dottheta_2*cos(theta_2))/2)))/2 + (m_V*(2*L2*sin(theta_2)*(dotx + L1*dottheta_1*sin(theta_1) + L2*dottheta_2*sin(theta_2)) + 2*L2*cos(theta_2)*(L1*dottheta_1*cos(theta_1) + L2*dottheta_2*cos(theta_2))))/2 , t));
% ddtdL_dotbeta = simplify(diff(I_V*dotbeta ,t));
% 
% eq_x = simplify(ddtdL_dotx - dL_x);
% eq_theta1 = simplify(ddtdL_dottheta1 - dL_theta1);
% eq_theta2 = simplify(ddtdL_dottheta2 - dL_theta2);
% eq_beta = simplify(ddtdL_dotbeta - dL_beta);

%% Pretty equations of motion

clear

syms m_V m_L1 m_L2 m_L3 I_V I_L2 I_L3 L1 L2 g_0
syms theta_1 theta_2 beta x dottheta_1 dottheta_2 dotbeta dotx ddottheta_1 ddottheta_2 ddotbeta ddotx
q = [theta_1 theta_2 beta x dottheta_1 dottheta_2 dotbeta dotx ddottheta_1 ddottheta_2 ddotbeta ddotx];

eq_x = m_L1*ddotx + m_L2*ddotx + m_L3*ddotx + m_V*ddotx + (L1*m_L2*sin(theta_1)*ddottheta_1)/2 + L1*m_L3*sin(theta_1)*ddottheta_1 + (L2*m_L3*sin(theta_2)*ddottheta_2)/2 + L1*m_V*sin(theta_1)*ddottheta_1 + L2*m_V*sin(theta_2)*ddottheta_2 + (L1*m_L2*cos(theta_1)*dottheta_1*dottheta_1)/2 + L1*m_L3*cos(theta_1)*dottheta_1*dottheta_1 + (L2*m_L3*cos(theta_2)*dottheta_2*dottheta_2)/2 + L1*m_V*cos(theta_1)*dottheta_1*dottheta_1 + L2*m_V*cos(theta_2)*dottheta_2*dottheta_2;
eq_theta1 = I_L2*ddottheta_1 + (L1*(g_0*m_L2*cos(theta_1) + 2*g_0*m_L3*cos(theta_1) + 2*g_0*m_V*cos(theta_1) - dotx*dottheta_1*m_L2*cos(theta_1) - 2*dotx*dottheta_1*m_L3*cos(theta_1) - 2*dotx*dottheta_1*m_V*cos(theta_1) + L2*dottheta_1*dottheta_2*m_L3*sin(theta_1 - theta_2) + 2*L2*dottheta_1*dottheta_2*m_V*sin(theta_1 - theta_2)))/2 + (L1^2*m_L2*ddottheta_1)/4 + L1^2*m_L3*ddottheta_1 + L1^2*m_V*ddottheta_1 + (L1*m_L2*sin(theta_1)*ddotx)/2 + L1*m_L3*sin(theta_1)*ddotx + L1*m_V*sin(theta_1)*ddotx + (L1*L2*m_L3*cos(theta_1 - theta_2)*ddottheta_2)/2 + L1*L2*m_V*cos(theta_1 - theta_2)*ddottheta_2 + (L1*m_L2*cos(theta_1)*dotx*dottheta_1)/2 + L1*m_L3*cos(theta_1)*dotx*dottheta_1 + L1*m_V*cos(theta_1)*dotx*dottheta_1 - (L1*L2*m_L3*sin(theta_1 - theta_2)*dottheta_2*dottheta_1)/2 + (L1*L2*m_L3*sin(theta_1 - theta_2)*dottheta_2*dottheta_2)/2 - L1*L2*m_V*sin(theta_1 - theta_2)*dottheta_2*dottheta_1 + L1*L2*m_V*sin(theta_1 - theta_2)*dottheta_2*dottheta_2;
eq_theta2 = I_L3*ddottheta_2 + (L2^2*m_L3*ddottheta_2)/4 + L2^2*m_V*ddottheta_2 - (L2*(m_L3 + 2*m_V)*(dotx*dottheta_2*cos(theta_2) - g_0*cos(theta_2) + L1*dottheta_1*dottheta_2*sin(theta_1 - theta_2)))/2 + (L2*m_L3*sin(theta_2)*ddotx)/2 + L2*m_V*sin(theta_2)*ddotx + (L1*L2*m_L3*cos(theta_1 - theta_2)*ddottheta_1)/2 + L1*L2*m_V*cos(theta_1 - theta_2)*ddottheta_1 + (L2*m_L3*cos(theta_2)*dotx*dottheta_2)/2 + L2*m_V*cos(theta_2)*dotx*dottheta_2 - (L1*L2*m_L3*sin(theta_1 - theta_2)*dottheta_1*dottheta_1)/2 + (L1*L2*m_L3*sin(theta_1 - theta_2)*dottheta_1*dottheta_2)/2 - L1*L2*m_V*sin(theta_1 - theta_2)*dottheta_1*dottheta_1 + L1*L2*m_V*sin(theta_1 - theta_2)*dottheta_1*dottheta_2;
eq_beta = I_V*ddotbeta;

eq_x = collect(eq_x, q);
eq_theta1 = collect(eq_theta1, q);
eq_theta2 = collect(eq_theta2, q);
eq_beta = collect(eq_beta, q);