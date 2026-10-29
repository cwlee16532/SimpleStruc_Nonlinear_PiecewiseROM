function [t, u_hist, maxTopDisp, regHist] = newmark_100dof_trilinear(M, C, params, dt, t, ag)
% newmark_100dof_trilinear : 100 DOF 전단구조물의
%                            3구간 선형 스토리 강성을 이용한
%                            Newmark-beta 지진해석 (누적 손상형)
%
%   [t, u_hist, maxTopDisp, regHist] = newmark_100dof_trilinear(M, C, params, dt, t, ag)
%
% 출력:
%   t         : 시간벡터 (1 x nt)
%   u_hist    : 각 층 변위 시계열 (N x nt)
%   maxTopDisp: 최상층 변위의 절대값 최대값
%   regHist   : 각 시간 step마다 각 스토리의 강성 구간 index (N x nt)
%               1 -> k1,  2 -> k2,  3 -> k3
%
% 비고:
%   비선형 강성행렬 K는 "현재 변위"가 아니라
%   "시간이력 중 최대 상대변위" umax 를 기반으로 조립됨.
%   ⇒ 한 번 항복하면(unloading 해도) 강성이 되돌아가지 않음.

% ---------------- 기본 설정 ----------------
N   = params.N;
uy1 = params.uy1;
uy2 = params.uy2;
k1  = params.k1;
k2  = params.k2;
k3  = params.k3;

nt = length(t);

% 초기조건
u = zeros(N,1);
v = zeros(N,1);
a = zeros(N,1);

u_hist  = zeros(N, nt);
u_hist(:,1) = u;

regHist = zeros(N, nt);  % 구간 이력 저장

% 스토리별 "최대 상대변위" 초기화
umax = zeros(N,1);       % 처음에는 모두 0

% 영향벡터 (전 층 동일 지반가속도)
r = ones(N,1);

% Newmark-beta 계수
beta  = 1/4;
gamma = 1/2;

% 초기 강성행렬 및 구간 index (umax = 0 이므로 모든 스토리 k1)
[K, reg] = assembleK_nonlinear(umax, N, uy1, uy2, k1, k2, k3);
regHist(:,1) = reg;

% 초기 가속도 (t = t(1) 에서 평형조건)
p0 = -M * r * ag(1);
a  = M \ (p0 - C*v - K*u);

% ---------------- 시간적분 루프 ----------------
for n = 1:nt-1
    
    % 현재 K는 umax (즉, 지금까지의 최대 상대변위)에 의해 결정된 상태
    
    % Newmark 예측 단계 (predictor)
    u_pred = u + dt*v + dt^2*(0.5 - beta)*a;
    v_pred = v + dt*(1 - gamma)*a;
    
    % 유효강성행렬
    Aeff = M + gamma*dt*C + beta*dt^2*K;
    
    % 하중 벡터 (지반가속도)
    p = -M * r * ag(n+1);
    
    % 잔차 (유효하중)
    rhs = p - C*v_pred - K*u_pred;
    
    % 새로운 가속도
    a_new = Aeff \ rhs;
    
    % 보정 단계 (corrector)
    u_new = u_pred + beta*dt^2*a_new;
    v_new = v_pred + gamma*dt*a_new;
    
    % 상태 업데이트
    u = u_new;
    v = v_new;
    a = a_new;
    
    % 변위 기록
    u_hist(:, n+1) = u;
    
    % ----- 여기서 umax 업데이트 (새로운 변위까지 반영) -----
    for i = 1:N
        if i == 1
            ue = u(1);          % 1층 - 기초(0층=0)
        else
            ue = u(i) - u(i-1); % i층 - (i-1)층
        end
        umax(i) = max(umax(i), abs(ue));  % 누적 최대 절대값
    end
    
    % 새롭게 업데이트된 umax 기반으로 강성행렬 및 구간 index 계산
    [K, reg] = assembleK_nonlinear(umax, N, uy1, uy2, k1, k2, k3);
    regHist(:, n+1) = reg;
    
end

maxTopDisp = max(abs(u_hist(N,:)));

end

