function [t_out, u_hist_rom, maxTop_rom, K_status_hist] = ...
    newmark_100dof_piecewise_eta_ROM( ...
        M, C, k_story, damageMask, eta1, eta2, dt, t_in, ag, th1, th2, Phi_rom)

% ROM 차수
r = size(Phi_rom,2);
N = size(M,1);
nt = length(t_in);

% ----- 전역 K (3상태) 구성 -----
k_story_vec_und = k_story * ones(N,1);
k_story_vec_d1  = k_story_vec_und;
k_story_vec_d2  = k_story_vec_und;

k_story_vec_d1(damageMask==1) = eta1 * k_story_vec_d1(damageMask==1);
k_story_vec_d2(damageMask==1) = eta2 * k_story_vec_d2(damageMask==1);

K_und = assembleK_from_story_stiffness(N, k_story_vec_und);
K_d1  = assembleK_from_story_stiffness(N, k_story_vec_d1);
K_d2  = assembleK_from_story_stiffness(N, k_story_vec_d2);

% ----- ROM 행렬 -----
Mr = Phi_rom' * M * Phi_rom;
Cr = Phi_rom' * C * Phi_rom;
Kr_und = Phi_rom' * K_und * Phi_rom;
Kr_d1  = Phi_rom' * K_d1  * Phi_rom;
Kr_d2  = Phi_rom' * K_d2  * Phi_rom;

% 기초가속도하중: p = -M*Gamma*ag(t)
Gamma = ones(N,1);
Fr = - (Phi_rom' * M * Gamma);   % scalar ag(t)에 대한 reduced force 형상 (r x 1)

% Newmark 계수
beta  = 1/4;
gamma = 1/2;

a0 = 1 / (beta * dt^2);
a1 = gamma / (beta * dt);
a2 = 1 / (beta * dt);
a3 = 1 / (2*beta) - 1;
a4 = gamma / beta - 1;
a5 = dt * (gamma/(2*beta) - 1);

% 초기 조건 (reduced)
qr    = zeros(r, nt);
qrd   = zeros(r, nt);
qrdd  = zeros(r, nt);

t_out = t_in;

% 초기 상태: undamaged
Kr_now = Kr_und;

% 초기 가속도
p0r = Fr * ag(1);  % (r x 1)
qrdd(:,1) = Mr \ (p0r - Cr*qrd(:,1) - Kr_now*qr(:,1));

% history
maxTop_run   = 0;
K_status_hist = zeros(1, nt);   % 1,2,3 상태 저장
K_status_hist(1) = 1;

for i = 1:nt-1
    
    % 현재까지의 roof disp (full로 복원해서 top DOF만)
    u_full_i = Phi_rom * qr(:,i);
    u_top_i  = u_full_i(N);
    maxTop_run = max(maxTop_run, abs(u_top_i));
    
    % 상태 결정
    if maxTop_run < th1
        Kr_now = Kr_und;
        status = 1;
    elseif maxTop_run < th2
        Kr_now = Kr_d1;
        status = 2;
    else
        Kr_now = Kr_d2;
        status = 3;
    end
    K_status_hist(i) = status;
    
    % 유효강성
    K_eff = Kr_now + a0*Mr + a1*Cr;
    
    % 외력 (i+1 시점)
    p_ext_r = Fr * ag(i+1);  % (r x 1)
    
    % 유효하중
    p_eff_r = p_ext_r ...
        + Mr*( a0*qr(:,i) + a2*qrd(:,i) + a3*qrdd(:,i) ) ...
        + Cr*( a1*qr(:,i) + a4*qrd(:,i) + a5*qrdd(:,i) );
    
    % 변위 (reduced)
    qr(:,i+1) = K_eff \ p_eff_r;
    
    % 가속도, 속도
    qrdd(:,i+1) = a0*(qr(:,i+1) - qr(:,i)) ...
                - a2*qrd(:,i) - a3*qrdd(:,i);
    qrd(:,i+1)  = qrd(:,i) + dt*((1-gamma)*qrdd(:,i) + gamma*qrdd(:,i+1));
end

% 마지막 스텝 상태 기록
K_status_hist(end) = K_status_hist(end-1);

% full DOF로 복원
u_hist_rom = Phi_rom * qr;

% 최대 top 변위
maxTop_rom = max(abs(u_hist_rom(N,:)));

end
