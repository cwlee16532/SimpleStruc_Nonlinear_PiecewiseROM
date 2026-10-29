function [t_out, u_hist, maxTopDisp, K_status_hist] = newmark_100dof_piecewise_eta( ...
    M, C, k_story, damageMask, eta1, eta2, dt, t, ag, th1, th2)
% newmark_100dof_piecewise_eta :
%   100 DOF 전단 모델에 대해 Newmark-beta(β=1/4, γ=1/2) 시간을력 해석 수행.
%   top 변위의 누적 최대값에 따라 1~20층 강성을 영구적으로 줄이는
%   piecewise-linear 파라미터화 모델.
%
%   상태 전환 규칙:
%     maxTop_run <  th1       -> undamaged   (eta = 1.0)
%     th1 <= maxTop_run < th2 -> damage 1    (eta = eta1)
%     maxTop_run >= th2       -> damage 2    (eta = eta2)
%
% 입력:
%   M, C       : 질량/감쇠행렬 (N x N)
%   k_story    : 기본 스토리강성 (undamaged)
%   damageMask : N x 1, 손상층 표시 (예: 1~20층 = 1, 나머지 0)
%   eta1,eta2  : damage state1, state2 에서의 강성 감소 계수
%   dt, t      : 시간간격, 시간벡터
%   ag         : 지반가속도 (1 x nt, [m/s^2])
%   th1, th2   : top 변위 threshold (예: 0.0065, 0.01)
%
% 출력:
%   t_out         : 시간벡터
%   u_hist        : 층 변위 시계열 (N x nt)
%   maxTopDisp    : 최상층 변위의 절대값 최대값
%   K_status_hist : 각 시간 step에서 상태 (1:undamaged, 2:damage1, 3:damage2)

N  = size(M,1);
nt = length(t);

% 초기조건
u = zeros(N,1);
v = zeros(N,1);
a = zeros(N,1);

u_hist = zeros(N, nt);
u_hist(:,1) = u;

K_status_hist = zeros(1, nt);   % 1,2,3 상태 기록
K_status = 1;                   % 초기 상태: undamaged
K_status_hist(1) = K_status;

% 영향벡터
r = ones(N,1);

% Newmark 계수
beta  = 1/4;
gamma = 1/2;

% 누적 top 변위 최대값
maxTop_run = 0;

% 초기 eta 선택
eta_use = 1.0;   % undamaged

% 초기 강성행렬
K = assembleK_linear_damaged(N, k_story, damageMask, eta_use);

% 초기 가속도
p0 = -M * r * ag(1);
a  = M \ (p0 - C*v - K*u);

for n = 1:nt-1
    
    % 현재까지의 max top 변위 업데이트
    maxTop_run = max(maxTop_run, abs(u(N)));
    
    % 상태 결정 (1/2/3)
    if maxTop_run < th1
        K_status = 1;   % undamaged
        eta_use  = 1.0;
    elseif maxTop_run < th2
        K_status = 2;   % damage 1
        eta_use  = eta1;
    else
        K_status = 3;   % damage 2
        eta_use  = eta2;
    end
    K_status_hist(n) = K_status;
    
    % 선택된 eta_use 에 따라 강성 재조립
    K = assembleK_linear_damaged(N, k_story, damageMask, eta_use);
    
    % 유효강성
    Aeff = M + gamma*dt*C + beta*dt^2*K;
    
    % predictor
    u_pred = u + dt*v + dt^2*(0.5 - beta)*a;
    v_pred = v + dt*(1 - gamma)*a;
    
    % 외력
    p = -M * r * ag(n+1);
    
    % 유효하중
    rhs = p - C*v_pred - K*u_pred;
    
    % 새로운 가속도
    a_new = Aeff \ rhs;
    
    % corrector
    u_new = u_pred + beta*dt^2*a_new;
    v_new = v_pred + gamma*dt*a_new;
    
    % 상태 업데이트
    u = u_new;
    v = v_new;
    a = a_new;
    
    % 기록
    u_hist(:, n+1) = u;
end

% 마지막 step 상태 기록
K_status_hist(nt) = K_status;

t_out = t;
maxTopDisp = max(abs(u_hist(N,:)));

end
