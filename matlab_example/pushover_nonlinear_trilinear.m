function [u_top_po, F_top_po, regHist_po] = pushover_nonlinear_trilinear(N, params, Ftop_max, nStep_po, Pshape)
 % pushover_nonlinear_trilinear :
%   3구간 누적손상(tri-linear) 모델을 사용하는 비선형 정적 푸쉬오버 해석
%   - top DOF에 수평하중을 0 -> Ftop_max 까지 점진적으로 증가
%   - 각 단계에서 Newton 반복으로 평형 맞추고,
%     umax 기반으로 강성(손상)을 업데이트
%
% 출력:
%   u_top_po  : 각 단계별 roof/top 변위 (1 x nStep_po)
%   F_top_po  : 각 단계별 top 하중 (1 x nStep_po)
%   regHist_po: 각 단계별 스토리 강성 구간 index (N x nStep_po)
%               1 -> k1, 2 -> k2, 3 -> k3
Fshape = Pshape(:);  % N x 1

uy1 = params.uy1;
uy2 = params.uy2;
k1  = params.k1;
k2  = params.k2;
k3  = params.k3;

% 초기 상태
u    = zeros(N,1);        % 층 변위
umax = zeros(N,1);        % 스토리별 최대 상대변위(누적)

u_top_po   = zeros(1, nStep_po);
F_top_po   = zeros(1, nStep_po);
regHist_po = zeros(N, nStep_po);

% top DOF 하중 벡터 베이스
eN = zeros(N,1);
eN(N) = 1;

% 하중 증분
F_top_steps = linspace(0, Ftop_max, nStep_po);

maxIter = 20;
tol     = 1e-6;

for istep = 1:nStep_po
    
    F_top = F_top_steps(istep);
    F = F_top * Fshape;   % 전 자유도 하중벡터
    
    % Newton 반복
    for iter = 1:maxIter
        
        % 현재 umax에 따른 재료강성행렬
        [K, ~] = assembleK_nonlinear(umax, N, uy1, uy2, k1, k2, k3);
        
        % 잔차 R = F - K u
        R  = F - K*u;
        
        % 변위 증분
        du = K \ R;
        u_new = u + du;
        
        % 새 변위에 따른 umax 업데이트(임시)
        umax_new = umax;
        for i = 1:N
            if i == 1
                ue = u_new(1);          % 1층-기초
            else
                ue = u_new(i) - u_new(i-1);
            end
            umax_new(i) = max(umax_new(i), abs(ue));
        end
        
        % 수렴 체크
        if norm(du) / max(1e-9, norm(u_new)) < tol
            u    = u_new;
            umax = umax_new;
            break;
        end
        
        % 업데이트 계속
        u    = u_new;
        umax = umax_new;
    end
    
    % 수렴 후 상태 기록
    u_top_po(istep) = u(N);
    F_top_po(istep) = F_top;
    
    [~, reg] = assembleK_nonlinear(umax, N, uy1, uy2, k1, k2, k3);
    regHist_po(:, istep) = reg;
end

end

