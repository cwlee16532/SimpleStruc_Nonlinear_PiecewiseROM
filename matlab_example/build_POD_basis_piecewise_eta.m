function Phi_pod = build_POD_basis_piecewise_eta(M, C, k_story, damageMask, eta1, eta2)
% build_POD_basis_piecewise_eta
% 100DOF 전단 구조물 + piecewise-eta 손상모델에 대해
% FRF 기반 POD 기저(10모드)를 생성

N = size(M,1);
r = 10;    % POD 차수

% --- 스토리별 강성 벡터 생성 (baseline)
k_story_vec_und = k_story * ones(N,1);

% 손상구역(damageMask==1)에 대해 eta1, eta2 적용
k_story_vec_d1  = k_story_vec_und;
k_story_vec_d2  = k_story_vec_und;
k_story_vec_d1(damageMask==1) = eta1 * k_story_vec_d1(damageMask==1);
k_story_vec_d2(damageMask==1) = eta2 * k_story_vec_d2(damageMask==1);

% 전역 K 조립
K_und = assembleK_from_story_stiffness(N, k_story_vec_und);
K_d1  = assembleK_from_story_stiffness(N, k_story_vec_d1);
K_d2  = assembleK_from_story_stiffness(N, k_story_vec_d2);

% FRF 계산용 주파수 범위 (예: 0.1~20 Hz)
fVec = (0.1:0.1:20)';
wVec = 2*pi*fVec;
nW   = length(wVec);

% 기초가속도에 대한 유효하중 형태(여기서는 -M*Gamma를 한 방향 입력으로 가정)
Gamma = ones(N,1);
GF = M * Gamma;   % base input 형상 (크기는 중요하지 않음, 형상만 중요)

% 3개 상태(undamaged, d1, d2)의 FRF를 모두 모음
FRF_snapshots = [];

for state = 1:3
    switch state
        case 1
            K_now = K_und;
        case 2
            K_now = K_d1;
        case 3
            K_now = K_d2;
    end
    
    H_state = zeros(N, nW);
    for iw = 1:nW
        w = wVec(iw);
        % (K + i w C - w^2 M)^-1 * GF
        H_state(:,iw) = (K_now + 1i*w*C - (w^2)*M) \ GF;
    end
    
    FRF_snapshots = [FRF_snapshots, H_state]; %#ok<AGROW>
end

% 실수/허수 분리해서 정보 풍부하게 만들고 싶으면 다음처럼도 가능:
% FRF_snapshots = [real(FRF_snapshots), imag(FRF_snapshots)];

% POD (SVD)
[U_pod, ~, ~] = svds(FRF_snapshots, r, "largest");

Phi_pod = U_pod;   % N x r

end

% -----------------------------------------------
% helper: 스토리별 강성벡터로 전역 K 조립
% -----------------------------------------------
