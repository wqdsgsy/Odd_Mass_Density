mr_xy = ppo(:,4)./(-1e-6*(2*pi*ppo(:,1)).^2);
mr_yx = ppo(:,3)./(-1e-6*(2*pi*ppo(:,1)).^2);
mr_x = ppo(:,2)./(-1e-6*(2*pi*ppo(:,1)).^2)-2.19e-3;
mr_y = ppo(:,5)./(-1e-6*(2*pi*ppo(:,1)).^2)-2.19e-3;
Eff_I = ppo(:,6)./(-1e-3*(2*pi*ppo(:,1)).^2);
mr_I = Eff_I-1.9e-7;
%
E_b = 270e9;   %Young modulus of host beam
rho_b = 7850;   %Mass density of host beam
bb = 1e-3;     %thickness of host beam
hh = 20e-3;    %unit cell length of host beam
AA = bb*hh;    %cross-sectional area of host beam
DD = E_b*hh^3*bb/12;   %Bending stiffness of host beam
hh0 = 4.5e-3;
AA0 = bb*hh0;
DD_b0 = E_b*hh0^3*bb/12;
dd = 7.75e-3;
E_m = 270e9;
% E_m = 60.46e9;
DD_b2 = DD_b0*2;
E_ml = 60.46e9;
rho_bl = rho_b/hh*hh0*2;
% E_ml = 270e9;
% DD_b1 = 2*E_b*hh0^3*bb/12 + 2*dd^2*E_b*AA0;
% DD_b2 = 2*E_b*hh0^3*bb/12;
% DD_b0 = 89.80;
% DD_b1 = DD_b0;
% DD_b2 = DD_b0;
% DD_m = 89.80;
% DD_m0 = DD*0.114;
% DD_m = DD*0.83;
LL = 9e-3;
LL_l = 10e-3;
% LL_2 = 6.5e-3;
% LL_3 = 10e-3;
pp = 0;
% for omega = 23e3*2*pi
for omega = 3e3*2*pi:100*2*pi:25e3*2*pi
    pp = pp + 1;
    k_l = sqrt(rho_b*omega^2/E_b);
    k_b = (rho_b*AA*omega^2/DD)^(1/4);
    k_lm = sqrt(rho_b*omega^2/E_m);
    k_bm = (rho_b*AA0*omega^2/DD_b0)^(1/4); 
    k_lml = sqrt(rho_bl*omega^2/E_ml);
    m_x = mr_x(pp);
    m_y = mr_y(pp);
    m_xy = mr_xy(pp);
    m_yx = mr_yx(pp);
    m_I = -mr_I(pp);
%     m_x = 0;
%     m_y = 0;
%     m_xy = 0;
%     m_yx = 0;
%     m_I = 0;
    %
    xx = -LL;
    xx_l = -LL_l;
    M0_n = [exp(1i*k_l*xx_l) exp(-1i*k_l*xx_l) 0 0 0 0; .../
        1i*k_l*E_b*AA*exp(1i*k_l*xx_l) -1i*k_l*E_b*AA*exp(-1i*k_l*xx_l) .../
        0 0 0 0; .../
        0 0 exp(1i*k_b*xx) exp(-1i*k_b*xx) exp(k_b*xx) exp(-k_b*xx); .../
        0 0 1i*k_b*exp(1i*k_b*xx) -1i*k_b*exp(-1i*k_b*xx) .../
        k_b*exp(k_b*xx) -k_b*exp(-k_b*xx); .../
        0 0 -k_b^2*DD*exp(1i*k_b*xx) -k_b^2*DD*exp(-1i*k_b*xx) .../
        k_b^2*DD*exp(k_b*xx) k_b^2*DD*exp(-k_b*xx); .../
        0 0 -1i*k_b^3*DD*exp(1i*k_b*xx) 1i*k_b^3*DD*exp(-1i*k_b*xx) .../
        k_b^3*DD*exp(k_b*xx) -k_b^3*DD*exp(-k_b*xx)];
    %
    xx = LL;
    xx_l = LL_l;
    M0_p = [exp(1i*k_l*xx_l) exp(-1i*k_l*xx_l) 0 0 0 0; .../
        1i*k_l*E_b*AA*exp(1i*k_l*xx_l) -1i*k_l*E_b*AA*exp(-1i*k_l*xx_l) .../
        0 0 0 0; .../
        0 0 exp(1i*k_b*xx) exp(-1i*k_b*xx) exp(k_b*xx) exp(-k_b*xx); .../
        0 0 1i*k_b*exp(1i*k_b*xx) -1i*k_b*exp(-1i*k_b*xx) .../
        k_b*exp(k_b*xx) -k_b*exp(-k_b*xx); .../
        0 0 -k_b^2*DD*exp(1i*k_b*xx) -k_b^2*DD*exp(-1i*k_b*xx) .../
        k_b^2*DD*exp(k_b*xx) k_b^2*DD*exp(-k_b*xx); .../
        0 0 -1i*k_b^3*DD*exp(1i*k_b*xx) 1i*k_b^3*DD*exp(-1i*k_b*xx) .../
        k_b^3*DD*exp(k_b*xx) -k_b^3*DD*exp(-k_b*xx)];
    %
    xx = -LL;
    xx_l = -LL_l;
    M1_n = [exp(1i*k_lml*xx_l) exp(-1i*k_lml*xx_l) 0 0 0 0; .../
        1i*k_lml*E_ml*AA*exp(1i*k_lml*xx_l) -1i*k_lml*E_ml*AA*exp(-1i*k_lml*xx_l) .../
        0 0 0 0; .../
        0 0 exp(1i*k_bm*xx) exp(-1i*k_bm*xx) exp(k_bm*xx) exp(-k_bm*xx); .../
        0 0 1i*k_bm*exp(1i*k_bm*xx) -1i*k_bm*exp(-1i*k_bm*xx) .../
        k_bm*exp(k_bm*xx) -k_bm*exp(-k_bm*xx); .../
        0 0 -k_bm^2*DD_b2*exp(1i*k_bm*xx) -k_bm^2*DD_b2*exp(-1i*k_bm*xx) .../
        k_bm^2*DD_b2*exp(k_bm*xx) k_bm^2*DD_b2*exp(-k_bm*xx); .../
        0 0 -1i*k_bm^3*DD_b2*exp(1i*k_bm*xx) 1i*k_bm^3*DD_b2*exp(-1i*k_bm*xx) .../
        k_bm^3*DD_b2*exp(k_bm*xx) -k_bm^3*DD_b2*exp(-k_bm*xx)];
    N1_n = -[0 0 0 0 0 0; .../
        0 0 0 0 0 0; .../
        0 0 0 0 0 0; .../
        0 0 0 0 0 0; .../
        0 0 0 0 1i*E_m*AA0*dd*k_lm*exp(1i*k_lm*xx) -1i*E_m*AA0*dd*k_lm*exp(-1i*k_lm*xx); .../
        0 0 0 0 -E_m*AA0*dd*k_lm^2*exp(1i*k_lm*xx) -E_m*AA0*dd*k_lm^2*exp(-1i*k_lm*xx)];
%     xx = -LL_2;
%     xx_l = -LL_3;
%     Q1 = [exp(1i*k_lml*xx_l) 0 0 0 0 0; .../
%         0 exp(-1i*k_lml*xx_l) 0 0 0 0; .../
%         0 0 exp(1i*k_bm*xx) 0 0 0; .../
%         0 0 0 exp(-1i*k_bm*xx) 0 0; .../
%         0 0 0 0 exp(k_bm*xx) 0; .../
%         0 0 0 0 0 exp(-k_bm*xx)];
    %
    xx = 0;
    M1_0 = [exp(1i*k_lml*xx) exp(-1i*k_lml*xx) 0 0 0 0; .../
        1i*k_lml*E_ml*AA*exp(1i*k_lml*xx) -1i*k_lml*E_ml*AA*exp(-1i*k_lml*xx) .../
        0 0 0 0; .../
        0 0 exp(1i*k_bm*xx) exp(-1i*k_bm*xx) exp(k_bm*xx) exp(-k_bm*xx); .../
        0 0 1i*k_bm*exp(1i*k_bm*xx) -1i*k_bm*exp(-1i*k_bm*xx) .../
        k_bm*exp(k_bm*xx) -k_bm*exp(-k_bm*xx); .../
        0 0 -k_bm^2*DD_b2*exp(1i*k_bm*xx) -k_bm^2*DD_b2*exp(-1i*k_bm*xx) .../
        k_bm^2*DD_b2*exp(k_bm*xx) k_bm^2*DD_b2*exp(-k_bm*xx); .../
        0 0 -1i*k_bm^3*DD_b2*exp(1i*k_bm*xx) 1i*k_bm^3*DD_b2*exp(-1i*k_bm*xx) .../
        k_bm^3*DD_b2*exp(k_bm*xx) -k_bm^3*DD_b2*exp(-k_bm*xx)];
%     N1_0 = -[0 0 0 0 0 0; .../
%         0 0 0 0 0 0; .../
%         0 0 0 0 0 0; .../
%         0 0 0 0 0 0; .../
%         0 0 0 0 1i*E_m*AA0*dd*k_lm*exp(1i*k_lm*xx) -1i*E_m*AA0*dd*k_lm*exp(-1i*k_lm*xx); .../
%         0 0 0 0 -E_m*AA0*dd*k_lm^2*exp(1i*k_lm*xx) -E_m*AA0*dd*k_lm^2*exp(-1i*k_lm*xx)];
    N1_0 = -[0 0 0 0 0 0; .../
        0 0 0 0 0 0; .../
        0 0 0 0 0 0; .../
        0 0 0 0 0 0; .../
        0 0 0 0 0 0; .../
        0 0 0 0 -E_m*AA0*dd*k_lm^2*exp(1i*k_lm*xx) -E_m*AA0*dd*k_lm^2*exp(-1i*k_lm*xx)];
    %
    xx = LL;
    xx_l = LL_l;
    M1_p = [exp(1i*k_lml*xx_l) exp(-1i*k_lml*xx_l) 0 0 0 0; .../
        1i*k_lml*E_ml*AA*exp(1i*k_lml*xx_l) -1i*k_lml*E_ml*AA*exp(-1i*k_lml*xx_l) .../
        0 0 0 0; .../
        0 0 exp(1i*k_bm*xx) exp(-1i*k_bm*xx) exp(k_bm*xx) exp(-k_bm*xx); .../
        0 0 1i*k_bm*exp(1i*k_bm*xx) -1i*k_bm*exp(-1i*k_bm*xx) .../
        k_bm*exp(k_bm*xx) -k_bm*exp(-k_bm*xx); .../
        0 0 -k_bm^2*DD_b2*exp(1i*k_bm*xx) -k_bm^2*DD_b2*exp(-1i*k_bm*xx) .../
        k_bm^2*DD_b2*exp(k_bm*xx) k_bm^2*DD_b2*exp(-k_bm*xx); .../
        0 0 -1i*k_bm^3*DD_b2*exp(1i*k_bm*xx) 1i*k_bm^3*DD_b2*exp(-1i*k_bm*xx) .../
        k_bm^3*DD_b2*exp(k_bm*xx) -k_bm^3*DD_b2*exp(-k_bm*xx)];
    N1_p = -[0 0 0 0 0 0; .../
        0 0 0 0 0 0; .../
        0 0 0 0 0 0; .../
        0 0 0 0 0 0; .../
        0 0 0 0 1i*E_m*AA0*dd*k_lm*exp(1i*k_lm*xx) -1i*E_m*AA0*dd*k_lm*exp(-1i*k_lm*xx); .../
        0 0 0 0 -E_m*AA0*dd*k_lm^2*exp(1i*k_lm*xx) -E_m*AA0*dd*k_lm^2*exp(-1i*k_lm*xx)];
%     xx = LL_2;
%     xx_l = LL_3;
%     Q2 = [exp(1i*k_lml*xx_l) 0 0 0 0 0; .../
%         0 exp(-1i*k_lml*xx_l) 0 0 0 0; .../
%         0 0 exp(1i*k_bm*xx) 0 0 0; .../
%         0 0 0 exp(-1i*k_bm*xx) 0 0; .../
%         0 0 0 0 exp(k_bm*xx) 0; .../
%         0 0 0 0 0 exp(-k_bm*xx)];
    %
%     TT1 = [1 1 -1 -1; .../
%         2i*k_lm*E_m*AA0*dd^2-omega^2*m_I -2i*k_lm*E_m*AA0*dd^2-omega^2*m_I .../
%         -2i*k_lm*E_m*AA0*dd^2 2i*k_lm*E_m*AA0*dd^2; .../
%         exp(-1i*k_lm*LL) exp(1i*k_lm*LL) 0 0; .../
%         0 0 exp(1i*k_lm*LL) exp(-1i*k_lm*LL)];
    TT1 = [1 1 -1 -1; .../
        1i*k_lm*E_m*AA0*dd^2 -1i*k_lm*E_m*AA0*dd^2 .../
        -1i*k_lm*E_m*AA0*dd^2 1i*k_lm*E_m*AA0*dd^2; .../
        exp(-1i*k_lm*LL) exp(1i*k_lm*LL) 0 0; .../
        0 0 exp(1i*k_lm*LL) exp(-1i*k_lm*LL)];
%     TT1 = [1 1 -1 -1; .../
%         1 1 0 0; .../
%         exp(-1i*k_lm*LL) exp(1i*k_lm*LL) 0 0; .../
%         0 0 exp(1i*k_lm*LL) exp(-1i*k_lm*LL)];
%     TT2 = dd*[0 0 0 0 0 0 0 0; .../
%         1i*k_bm*exp(-1i*k_bm*LL)*omega^2*m_I -1i*k_bm*exp(1i*k_bm*LL)*omega^2*m_I .../
%         k_bm*exp(-k_bm*LL)*omega^2*m_I -k_bm*exp(k_bm*LL)*omega^2*m_I .../
%         1i*k_bm*exp(1i*k_bm*LL)*omega^2*m_I -1i*k_bm*exp(-1i*k_bm*LL)*omega^2*m_I .../
%         k_bm*exp(k_bm*LL)*omega^2*m_I -k_bm*exp(-k_bm*LL)*omega^2*m_I; .../
%         -1i*k_bm*exp(-1i*k_bm*LL) 1i*k_bm*exp(1i*k_bm*LL) .../
%         -k_bm*exp(-k_bm*LL) k_bm*exp(k_bm*LL) 0 0 0 0; .../
%         0 0 0 0 -1i*k_bm*exp(1i*k_bm*LL) 1i*k_bm*exp(-1i*k_bm*LL) .../
%         -k_bm*exp(k_bm*LL) k_bm*exp(-k_bm*LL)];
    TT2 = dd*[0 0 0 0 0 0 0 0; .../
        exp(-1i*k_bm*LL)/LL*omega^2*m_I exp(1i*k_bm*LL)/LL*omega^2*m_I .../
        exp(-k_bm*LL)/LL*omega^2*m_I exp(k_bm*LL)/LL*omega^2*m_I .../
        -exp(1i*k_bm*LL)/LL*omega^2*m_I -exp(-1i*k_bm*LL)/LL*omega^2*m_I .../
        -exp(k_bm*LL)/LL*omega^2*m_I -exp(-k_bm*LL)/LL*omega^2*m_I; .../
        -1i*k_bm*exp(-1i*k_bm*LL) 1i*k_bm*exp(1i*k_bm*LL) .../
        -k_bm*exp(-k_bm*LL) k_bm*exp(k_bm*LL) 0 0 0 0; .../
        0 0 0 0 -1i*k_bm*exp(1i*k_bm*LL) 1i*k_bm*exp(-1i*k_bm*LL) .../
        -k_bm*exp(k_bm*LL) k_bm*exp(-k_bm*LL)];
    TTn = inv(TT1)*TT2;
    TTn_f11(1:6, 1:6) = 0;
    TTn_f11(5:6, 3:6) = TTn(1:2, 1:4);
    TTn_f12(1:6, 1:6) = 0;
    TTn_f12(5:6, 3:6) = TTn(1:2, 5:8);
    TTn_f21(1:6, 1:6) = 0;
    TTn_f21(5:6, 3:6) = TTn(3:4, 1:4);
    TTn_f22(1:6, 1:6) = 0;
    TTn_f22(5:6, 3:6) = TTn(3:4, 5:8);
    %
    FF = [0 0 0 0 0 0;
        omega^2*m_x omega^2*m_x omega^2*m_xy .../
        omega^2*m_xy omega^2*m_xy omega^2*m_xy; .../
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        -omega^2*m_yx -omega^2*m_yx -omega^2*m_y .../
        -omega^2*m_y -omega^2*m_y -omega^2*m_y];
    %
    M1_n = M1_n + 2*N1_n*TTn_f11;
    M1_p = M1_p + 2*N1_p*TTn_f22;
    M1_0n1 = M1_0 + 2*N1_0*TTn_f11 - 2*N1_0*TTn_f21 - FF;
    M1_0n2 = M1_0 + 2*N1_0*TTn_f22 - 2*N1_0*TTn_f12;
    %
%     QQ = Q2;
%     QQ = Q1 + Q2;
%     TT = inv(M0_n)*M1_n*inv(M1_0+FF)*M1_0*inv(M1_p)*M0_p;
    TT = inv(M0_n)*(M1_n*inv(M1_0n1)*M1_0n2+2*N1_n*TTn_f12)* .../
        inv(M1_p+2*N1_p*TTn_f21*inv(M1_0n1)*M1_0n2)*M0_p;
%     TT = inv(M0_n)*(M1_n+FF*QQ/4)*inv(M1_p-FF*QQ/4)*M0_p;
    M_l = [TT(1,1) TT(1,3) TT(1,6) 0 0 0; .../
        TT(3,1) TT(3,3) TT(3,6) 0 0 0; .../
        TT(6,1) TT(6,3) TT(6,6) 0 0 0; .../
        TT(2,1) TT(2,3) TT(2,6) -1 0 0; .../
        TT(4,1) TT(4,3) TT(4,6) 0 -1 0; .../
        TT(5,1) TT(5,3) TT(5,6) 0 0 -1];
    %
    V_r = [1; 0; 0; 0; 0; 0];
%     V_r = [0; 1; 0; 0; 0; 0];
    %
    T_t = M_l\V_r;
    freq_p(pp) = omega/2/pi;
    tt_l(pp) = T_t(1);
    tt_t(pp) = T_t(2);
end
plot(freq_p/1e3, abs(tt_l))
hold on
plot(freq_p/1e3, abs(tt_t))
hold on

