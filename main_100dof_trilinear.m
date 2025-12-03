%% main_100dof_trilinear.m
% 100 DOF 전단구조물 + 3구간 선형 스토리강성 + Newmark-beta 지진해석
% - 실제 지진기록 (경주/포항) CSV 파일 사용
clear; clc; close all;

%% 기본 설정
N = 100;           % 자유도(층수)
m_story = 10000;     % 층당 lumped mass (단위질량으로 설정)
M = m_story * eye(N);

%% 1차 고유진동수 ~ 4 Hz가 되도록 스토리 강성 튜닝 (k_story는 이미 맞춰둔 값 사용)
k_story = 2.5857e+10;   % 사용자가 조정해 둔 값

K0 = assembleK_linear(N, k_story);

% 일반화 고유치문제 풀어서 1, 2차 고유진동수 확인
[phi, D] = eig(K0, M);
w = sqrt(diag(D));
[~, idx] = sort(w);
w_sorted = w(idx);
w1 = w_sorted(1);
w2 = w_sorted(2);
f1 = w1 / (2*pi);
f2 = w2 / (2*pi);
f = w_sorted / (2*pi);

fprintf('1차 고유진동수 f1 = %.3f Hz\n', f1);
fprintf('2차 고유진동수 f2 = %.3f Hz\n', f2);

%% Rayleigh 감쇠행렬 (ξ = 5% at 1, 2차 모드)
zeta = 0.05;
A = [1/(2*w1), w1/2;
     1/(2*w2), w2/2];
x = A \ [zeta; zeta];
alpha = x(1);
beta  = x(2);

C = alpha*M + beta*K0;

%% 3구간 선형 스토리 모델 파라미터
% 변위 경계 (필요에 따라 수정 가능)
uy1 = 0.0001;   % 1차 경계 변위 (예: 0.5cm)
uy2 = 0.0002;    % 2차 경계 변위 (예: 2cm)

% 강성 (손상 전/후를 단순 모델링)
k1 = k_story;           % 탄성영역 강성
k2 = 0.5 * k_story;     % 1차 손상 후 강성
k3 = 0.1 * k_story;     % 심한 손상 후 강성

params.N   = N;
params.uy1 = uy1;
params.uy2 = uy2;
params.k1  = k1;
params.k2  = k2;
params.k3  = k3;

%% 3선형 스토리 강성 모델 backbone 곡선 플롯 (논문용 예시)

% 플롯할 최대 변위 (양쪽으로 uy2의 1.5배 정도까지)
u_max = 1.5 * uy2;
u = linspace(-u_max, u_max, 400);   % 대칭 backbone 그리기
F = zeros(size(u));

% tri-linear 힘–변위 관계 정의
% 구간 1: |u| <= uy1   → 기울기 k1
% 구간 2: uy1 < |u| <= uy2 → 기울기 k2 (감소강성)
% 구간 3: |u| >  uy2       → 기울기 k3 (심한 손상)

idx1 = abs(u) <= uy1;
idx2 = abs(u) > uy1 & abs(u) <= uy2;
idx3 = abs(u) > uy2;

% 구간 1
F(idx1) = k1 * u(idx1);

% 구간 2 (연속성 보장: 첫 구간의 끝값 + k2*(Δu))
F_y1 = k1 * uy1;  % 첫 항복점 힘
F(idx2) = sign(u(idx2)) .* ( F_y1 + k2 * (abs(u(idx2)) - uy1) );

% 구간 3 (마찬가지로 연속)
F_y2 = F_y1 + k2 * (uy2 - uy1);  % 두 번째 경계점 힘
F(idx3) = sign(u(idx3)) .* ( F_y2 + k3 * (abs(u(idx3)) - uy2) );

% 플롯
figure;
plot(u*1000, F/1e6, 'LineWidth', 1.8);  % [mm, MN] 스케일 예시
hold on;

% 양의 방향 경계점 마커 표시
plot([-uy1 -uy2 uy1 uy2]*1000, [-F_y1 -F_y2 F_y1 F_y2]/1e6, 'ko', 'MarkerFaceColor', 'k');

% 축/레이블/스타일 (논문용)
xlabel('Story displacement, u (mm)');
ylabel('Story shear force, F (MN)');
title('Tri-linear story stiffness model');
grid on;
box on;

% 축 범위 약간 여유
xlim([-u_max u_max]*1000);

% 글꼴/선두께 등 논문용 옵션 (필요시)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1);

legend({'Tri-linear backbone', 'Yield points (u_{y1}, u_{y2})'}, ...
       'Location', 'NorthWest');
set(gcf, 'Color', 'w');       % figure 배경 흰색
set(gcf,'Units','inches','Position',[1 1 6 4])  % 6x4 inch figure


%% 지진기록 불러오기 (경주 또는 포항 선택)
% eqFile = 'GY_USN_3axes.csv';   % 경주지진
eqFile = 'Po_PHA2_3axes.csv';    % 포항지진

eqTable = readtable(eqFile);
eq      = eqTable{:,:};        % numeric 배열로 변환

time_raw = eq(:,1);            % [s]
acc_H1   = eq(:,2);            % 첫 번째 수평성분 (H1)

% baseline offset 제거
acc_H1 = acc_H1 - acc_H1(1);

% 시간 정보
dt   = time_raw(2) - time_raw(1);
tEnd = time_raw(end) - time_raw(1);
t    = (time_raw - time_raw(1))';   % 1 x nt 로 사용 (row vector)
nt   = numel(t);

% 원기록 지반가속도 (단위: m/s^2 라고 가정)
ag0 = acc_H1';       % 1 x nt

% 원기록 PGA
g = 9.81;
PGA0 = max(abs(ag0));          % [m/s^2]
PGA0_g = PGA0 / g;             % [g]

fprintf('원기록 PGA = %.3f g\n', PGA0_g);

%% 원 지반가속도 시간이력 플롯 (논문용)

figure;
plot(t, ag0/g, 'k-', 'LineWidth', 1.2);  % [g] 단위로 플롯
grid on; box on;

xlabel('Time (s)');
ylabel('Ground acceleration (g)');
% title(sprintf('Recorded ground motion (station: %s)', eqFile), ...
      % 'Interpreter','none');

% 축 범위 및 표시 옵션 (필요시 조정)
xlim([t(1) t(end)]);
% ylim([-1.1 1.1]*max(abs(ag0/g)));  % 여유 있게 설정하고 싶으면 사용

set(gca, 'FontName','Times New Roman', 'FontSize', 12, 'LineWidth', 1);

% PGA 값 표시용 텍스트 (옵션)
text(0.05*(t(end)-t(1)) + t(1), -0.95*PGA0_g, ...
     sprintf('PGA = %.3f g', PGA0_g), ...
     'FontName','Times New Roman', 'FontSize', 11);
set(gcf, 'Color', 'w');       % figure 배경 흰색
set(gcf,'Units','inches','Position',[1 1 6 4])  % 6x4 inch figure

%% 스케일링 리스트 (예: 원기록의 0.2 ~ 1.0배)
sf_list = 0:0.025:0.5;         % 필요에 따라 변경
ncase   = numel(sf_list);
maxTop  = zeros(ncase,1);
PGA_list = sf_list;   % 각 케이스 PGA (g 단위)

% 마지막 케이스(time history & 항복이력 plotting 위해 저장)
u_hist_last  = [];
t_out_last   = [];
regHist_last = [];
regHist_data = [];
for ic = 1:ncase
    sf = sf_list(ic);
    ag = sf * ag0/PGA0*g;   % 스케일된 지반가속도 [m/s^2]
    
    [t_out, u_hist, maxTop(ic), regHist] = newmark_100dof_trilinear(M, C, params, dt, t, ag);
    
    fprintf('Case %2d: scale = %.2f,  PGA = %.3f g,  max top disp = %.5f m\n', ...
        ic, sf, PGA_list(ic), maxTop(ic));
    
    if ic == ncase
        u_hist_last  = u_hist;
        t_out_last   = t_out;
        regHist_last = regHist;
    end
    regHist_data{ic} = regHist;
    u_hist_data{ic} = u_hist;
end

%% 성능곡선 (PGA vs 최대 top disp)
figure;
plot(maxTop, sf_list,'-o', 'LineWidth', 1.5);
ylabel('PGA (g)');
xlabel('Max top displacement (m)');
title(sprintf('Performance curve (record: %s)', eqFile), 'Interpreter', 'none');
grid on;

%% 가장 큰 스케일 케이스(top 층 변위 시계열)
if ~isempty(u_hist_last)
    figure;
    % plot(t_out_last, u_hist_last(end,:), 'LineWidth', 1.2);
    plot(t_out_last, u_hist_data{15}(100,:), 'LineWidth', 1.2);
    xlabel('Time (s)');
    ylabel('Top floor displacement (m)');
    title(sprintf('Top displacement time history (record: %s, scale = %.2f, PGA = %.3f g)', ...
        eqFile, sf_list(end), PGA_list(end)), ...
        'Interpreter', 'none');
    grid on;
end

%% 요소별(스토리별) 강성 구간 이력 플롯 (1:k1, 2:k2, 3:k3)
if ~isempty(regHist_last)
    figure;
    imagesc(t_out_last, 1:N, regHist_data{19});  %11:0.250g, 15:0.350g, 19: 0.450g
    axis xy;
    xlabel('Time (s)');
    ylabel('Story / Element index');
    % title(sprintf('Yielding history (record: %s, scale = %.2f)', ...
    %     eqFile, sf_list(end)), ...
        % 'Interpreter', 'none');

    cb = colorbar;
    cb.Ticks = [1 2 3];
    cb.TickLabels = {'k1','k2','k3'};

    caxis([1 3]);

    % ===========================
    % custom colormap (blue / orange / red)
    % ===========================
    cmap = [
        0     0.4470 0.7410   % blue
        0.8500 0.3250 0.0980   % orange
        0.6350 0.0780 0.1840   % red
    ];
    colormap(cmap);
end
set(gcf, 'Color', 'w');       % figure 배경 흰색
set(gcf,'Units','inches','Position',[1 1 6 4])  % 6x4 inch figure
set(gca, 'FontName','Times New Roman', 'FontSize', 16, 'LineWidth', 1);

% %% ---------------------------------------------------------
% %% 비선형 푸쉬오버 해석 (top DOF 점진적 가력, 3구간 누적손상 모델 사용)
% %% ---------------------------------------------------------
% 
% % 1) 우선 초기 선형 K0로 단위 top 하중에 대한 변위 형상 계산
% eN = zeros(N,1);
% eN(N) = 1;                % top DOF에 1.0의 수평하중
% 
% u_unit = K0 \ eN;         % K0 * u = eN
% roof_disp_per_unit = u_unit(N);      % 단위 하중당 roof 변위
% k_roof = 1 / roof_disp_per_unit;     % 초기 등가 roof 강성
% 
% % 2) 목표 roof 변위 설정 (예: uy2의 3배 정도)
% u_target = 2 * uy2;
% 
% % 해당 roof 변위를 만들기 위한 top 하중 (초기 탄성 기준)
% Ftop_max = k_roof * u_target;
% nStep_po = 50;
% 
% % --- main에서 1차 모드 벡터를 사용해 가력 패턴 정의 ---
% [phi,~] = eig(K0,M);
% [~,idx] = sort(w);
% phi1 = phi(:, idx(1));         % 1차 모드
% Pshape = phi1 / max(phi1);     % 최상층 = 1 로 정규화
% 
% % pushover 호출을 이렇게 변경:
% [u_top_po, Fscale_po, regHist_po] = pushover_nonlinear_trilinear(N, params, Ftop_max, nStep_po, Pshape);
% 
% 
% % 4) 등가 평균 가속도(PGA 형태)로 변환 (일반 분포하중 버전)
% M_tot = sum(diag(M));      % 총 질량
% g     = 9.81;
% 
% % Pshape: N x 1 모드형/분포형 가력 패턴 (최상층 = 1 등)
% F_total_po = Fscale_po * sum(Pshape);   % 각 단계별 총 수평력 (base shear)
% a_eq   = F_total_po / M_tot;            % 평균 가속도 [m/s^2]
% PGA_eq = a_eq / g;                      % [g]
% 
% % 5) 비선형 푸쉬오버 성능곡선 (F_total - roof disp)
% figure;
% plot(u_top_po, F_total_po, '-o', 'LineWidth', 1.5);
% xlabel('Roof displacement (m)');
% ylabel('Total lateral force / base shear (N)');
% title('Nonlinear static pushover curve (V_b - roof disp)');
% grid on;
% 
% % 6) 동적 성능곡선과 비교하기 쉬운 형태 (Roof disp vs 등가 PGA)
% figure;
% plot(u_top_po, PGA_eq, '-o', 'LineWidth', 1.5); hold on;
% plot(maxTop, PGA_list, 's-', 'LineWidth', 1.5);
% xlabel('Roof / top displacement (m)');
% ylabel('PGA (g)');
% legend('Nonlinear static pushover (equiv. PGA)', ...
%        'Dynamic TH results', ...
%        'Location','NorthWest');
% title('Comparison of nonlinear pushover and dynamic performance');
% grid on;
% 
% % 7) 푸쉬오버 중 층별 강성구간(reg) 이력도 보고 싶다면:
% figure;
% imagesc(1:nStep_po, 1:N, regHist_po);   % x: pushover step index
% axis xy;
% xlabel('Pushover load step');
% ylabel('Story / Element index');
% title('Yielding history during nonlinear static pushover');
% cb = colorbar;
% cb.Ticks = [1 2 3];
% cb.TickLabels = {'k1','k2','k3'};
% caxis([1 3]);
% colormap(jet(3));

%% ---------------------------------------------------------
%% 1) 기준 tri-linear 모델에서 특정 PGA의 target top 변위 추출
%% ---------------------------------------------------------

% 기준 tri-linear 모델 결과:
%  - maxTop(ic) : 각 PGA에서 최상층 max disp (기존에서 이미 계산됨)
%  - PGA_list(ic): 각 케이스별 PGA [g] (sf_list와 동일하게 사용 중)

target_PGA1 = 0.350;   % [g]
target_PGA2 = 0.450;    % [g]

idx1 = find(abs(PGA_list - target_PGA1) < 1e-8, 1);
idx2 = find(abs(PGA_list - target_PGA2) < 1e-8, 1);

if isempty(idx1) || isempty(idx2)
    error('PGA_list 에 0.350 g 또는 0.45 g 가 포함되어 있지 않습니다. sf_list 범위를 확인하세요.');
end

u_ref_0350 = maxTop(idx1);   % baseline tri-linear 모델의 top max disp at 0.375 g
u_ref_0450  = maxTop(idx2);   % baseline tri-linear 모델의 top max disp at 0.45 g

fprintf('Baseline tri-linear: PGA=0.350g -> maxTop = %.5e m\n', u_ref_0350);
fprintf('Baseline tri-linear: PGA=0.450g  -> maxTop = %.5e m\n', u_ref_0450);

%% ---------------------------------------------------------
%% 2) piecewise-eta 모델에서 eta1, eta2 캘리브레이션
%%    - 상태 전환 규칙: maxTop_run 기준
%%      maxTop < 0.0065  -> eta = 1.0 (undamaged)
%%      0.0065 <= maxTop < 0.01 -> eta = eta1 (damage1)
%%      maxTop >= 0.01          -> eta = eta2 (damage2)
%% ---------------------------------------------------------

th1 = 0.0065;    % 첫 번째 damage threshold
th2 = 0.009;    % 두 번째 damage threshold

% 1~20층을 damage zone으로 설정
damageMask = [ones(20,1); zeros(N-20,1)];

% eta 후보 grid
eta_grid = linspace(0.02, 0.5, 500);

% ---- (a) eta1 캘리브레이션: PGA = 0.375g 에서 top disp = u_ref_0375 맞추기 ----
err1 = zeros(size(eta_grid));

for k = 1:length(eta_grid)
    eta_candidate = eta_grid(k);
    
    % 이 단계에서는 eta1=eta_candidate, eta2는 임시값(같거나 1.0)으로 두되
    % 0.375g에서 maxTop이 주로 th1~th2 사이에 있을 것이므로 eta2 영향은 거의 없음
    eta1_test = eta_candidate;
    eta2_test = eta_candidate;  % 혹은 1.0 등으로 두어도 무방
    
    sf = target_PGA1;
    ag_eta = sf * ag0 / PGA0 * g;   % 원기록 PGA를 sf[g]로 맞추는 scaling
    
    [~, u_hist_eta, maxTop_eta, ~] = newmark_100dof_piecewise_eta( ...
        M, C, k_story, damageMask, eta1_test, eta2_test, dt, t, ag_eta, th1, th2);
    
    err1(k) = abs(maxTop_eta - u_ref_0350);
end

[~, idx_eta1] = min(err1);
eta1 = eta_grid(idx_eta1);
fprintf('Calibrated eta1 (for 0.375g) = %.3f\n', eta1);

% ---- (b) eta2 캘리브레이션: PGA = 0.45g 에서 top disp = u_ref_045 맞추기 ----
err2 = zeros(size(eta_grid));

for k = 1:length(eta_grid)
    eta_candidate = eta_grid(k);
    
    eta1_test = eta1;            % 이제 eta1은 고정
    eta2_test = eta_candidate;   % eta2만 바꿔본다
    
    sf = target_PGA2;
    ag_eta = sf * ag0 / PGA0 * g;
    
    [~, u_hist_eta, maxTop_eta, ~] = newmark_100dof_piecewise_eta( ...
        M, C, k_story, damageMask, eta1_test, eta2_test, dt, t, ag_eta, th1, th2);
    
    err2(k) = abs(maxTop_eta - u_ref_0450);
end

[~, idx_eta2] = min(err2);
eta2 = eta_grid(idx_eta2);
fprintf('Calibrated eta2 (for 0.45g) = %.3f\n', eta2);

%% ---------------------------------------------------------
%% 3) 캘리브레이션된 piecewise-eta 모델로 전체 PGA 범위 성능곡선 계산
%% ---------------------------------------------------------

maxTop_param_pw   = zeros(size(maxTop));
K_status_all      = cell(ncase,1);  % 필요시 상태 이력 저장
u_hist_pw_data    = cell(ncase,1);  % piecewise-eta 모델 변위 시계열 저장

for ic = 1:ncase
    sf = sf_list(ic);   % 이게 곧 해당 케이스의 PGA [g]
    
    ag_eta = sf * ag0 / PGA0 * g;
    
    [t_pw, u_hist_pw, maxTop_pw, K_status_hist_pw] = newmark_100dof_piecewise_eta( ...
        M, C, k_story, damageMask, eta1, eta2, dt, t, ag_eta, th1, th2);
    
    maxTop_param_pw(ic) = maxTop_pw;
    K_status_all{ic}    = K_status_hist_pw;   % 필요 시 시각화 용도
    u_hist_pw_data{ic}  = u_hist_pw;         % time history 저장
    K_status_hist_pw_data{ic}  = K_status_hist_pw;         % time history 저장

end

%% ---------------------------------------------------------
%% 4) 기존 tri-linear 성능곡선 vs piecewise-eta 모델 성능곡선 비교 플롯
%% ---------------------------------------------------------

figure;
plot(maxTop,          PGA_list, 'bo-', 'LineWidth', 1.5); hold on;
plot(maxTop_param_pw, PGA_list, 'rs--', 'LineWidth', 1.5);
ylabel('PGA (g)');
xlabel('Max top displacement (m)');
legend('Baseline tri-linear model', 'Piecewise-eta parametrized model', ...
       'Location','NorthWest');
title(sprintf('Baseline vs piecewise-eta model performance (record: %s)', eqFile), ...
      'Interpreter', 'none');
grid on;
ylim([0 0.45])

%% ---------------------------------------------------------
%% 6) 캘리브레이션된 piecewise-eta 모델의 상태이력(K_status) 플롯
%%    (1: undamaged, 2: damage1, 3: damage2)
%% ---------------------------------------------------------

case_list  = [idx1, idx2];
case_PGAs  = [target_PGA1, target_PGA2];
%11:0.250g, 15:0.350g, 19: 0.450g
for ic = [15 19]

    K_status_hist = K_status_hist_pw_data{ic};  % 1 x nt

    figure;
    stairs(t, K_status_hist, 'LineWidth', 1.5);   % 계단형이 직관적
    xlabel('Time (s)');
    ylabel('System State');
    % title(sprintf('Piecewise-eta damage state history at PGA = %.3f g (case %d)', ...
    %     case_PGAs(jj), ic));

    yticks([1 2 3]);
yticklabels({'undamaged', ...
             'damaged #1', ...
             'damaged #2'});


    ylim([0.9 3.1]);
    grid on;
    set(gcf, 'Color', 'w');       % figure 배경 흰색
    set(gcf,'Units','inches','Position',[1 1 6 4])  % 6x4 inch figure
    set(gca, 'FontName','Times New Roman', 'FontSize', 12, 'LineWidth', 1);
end


%% ---------------------------------------------------------
%% 5) 같은 PGA 조건에서 top 변위 time history 비교 (tri-linear vs piecewise-eta)
%% ---------------------------------------------------------

% 0.350 g, 0.425 g 케이스 index는 위에서 이미 idx1, idx2 로 구해둔 상태
case_list  = [idx1, idx2];
case_PGAs  = [target_PGA1, target_PGA2];

for jj = 1:length(case_list)
    ic = case_list(jj);
    
    % tri-linear 모델의 top 변위
    u_top_tri = u_hist_data{ic}(N,:);          % N=100번째 DOF (top floor)
    
    % piecewise-eta 모델의 top 변위
    u_top_eta = u_hist_pw_data{ic}(N,:);       % 같은 시간축 t 사용
    
    figure;
    plot(t, u_top_tri, 'b-', 'LineWidth', 1.2); hold on;
    plot(t, u_top_eta, 'r--', 'LineWidth', 1.2);
    xlabel('Time (s)');
    ylabel('Top floor displacement (m)');
    title(sprintf('Top displacement time history at PGA = %.3f g (case %d)', ...
        case_PGAs(jj), ic));
    legend('Baseline tri-linear', 'Piecewise-eta model', ...
           'Location','NorthEast');
    grid on;
end

%% ----------------- FRF 기반 POD-ROM (10 modes) 생성 -----------------
Phi_pod = build_POD_basis_piecewise_eta(M, C, k_story, damageMask, eta1, eta2);

%% 3') POD-ROM 기반 piecewise-eta 모델로 전체 PGA 범위 성능곡선 계산

maxTop_param_pw_rom   = zeros(size(maxTop));
u_hist_pw_rom_data    = cell(ncase,1);

for ic = 1:ncase
    sf = sf_list(ic);   % 이게 곧 해당 케이스의 PGA [g]
    
    ag_eta = sf * ag0 / PGA0 * g;
    
    [t_pw_rom, u_hist_pw_rom, maxTop_pw_rom, K_status_hist_pw_rom] = ...
        newmark_100dof_piecewise_eta_ROM( ...
            M, C, k_story, damageMask, eta1, eta2, dt, t, ag_eta, th1, th2, Phi_pod);
    
    maxTop_param_pw_rom(ic) = maxTop_pw_rom;
    u_hist_pw_rom_data{ic}  = u_hist_pw_rom;   % time history 저장
end

%% ---------------------------------------------------------
%% 4) tri-linear vs piecewise-eta vs POD-ROM 성능곡선 비교
figure;
hold on;

% Baseline tri-linear
plot(maxTop, PGA_list, ...
    'o-', 'Color', 'k', 'LineWidth', 1.8, 'MarkerFaceColor','k');

% Piecewise-eta (FOM)
plot(maxTop_param_pw, PGA_list, ...
    's-', 'Color', [0 0.4470 0.7410], ...   % 진한 파랑
    'LineWidth', 1.8, 'MarkerFaceColor',[0 0.4470 0.7410]);

% POD-ROM
plot(maxTop_param_pw_rom, PGA_list, ...
    '^--' , 'Color', [0.8500 0.3250 0.0980], ...   % 주황-빨강
    'LineWidth', 1.8, 'MarkerFaceColor',[0.8500 0.3250 0.0980]);

ylabel('PGA (g)');
xlabel('Max top displacement (m)');
ylim([0 0.45])
legend('Baseline tri-linear model', ...
       'Piecewise linear parametrized model (FOM)', ...
       'POD-ROM piecewise linear model (10 modes)', ...
       'Location','best');

% title(sprintf('Baseline vs piecewise-eta vs POD-ROM (record: %s)', eqFile), ...
%       'Interpreter','none');

grid on;
box on;

set(gcf, 'Color', 'w');       % figure 배경 흰색
set(gcf,'Units','inches','Position',[1 1 6 4])  % 6x4 inch figure
set(gca, 'FontName','Times New Roman', 'FontSize', 12, 'LineWidth', 1);
%% ---------------------------------------------------------
%% 5) 같은 PGA 조건에서 top 변위 time history 비교 (3모델)
%% ---------------------------------------------------------

case_list  = [idx1, idx2];
case_PGAs  = [target_PGA1, target_PGA2];

for jj = 1:length(case_list)
    ic = case_list(jj);
    
    % tri-linear 모델
    u_top_tri = u_hist_data{ic}(N,:);           
    
    % piecewise-eta (FOM)
    u_top_eta = u_hist_pw_data{ic}(N,:);       
    
    % POD-ROM piecewise-eta
    u_top_rom = u_hist_pw_rom_data{ic}(N,:);   
    
    figure;
    plot(t, u_top_tri, 'b-',  'LineWidth', 1.2); hold on;
    plot(t, u_top_eta, 'r--', 'LineWidth', 1.2);
    plot(t, u_top_rom, 'k-.', 'LineWidth', 1.2);
    xlabel('Time (s)');
    ylabel('Top floor displacement (m)');
    title(sprintf('Top displacement time history at PGA = %.3f g (case %d)', ...
        case_PGAs(jj), ic));
    legend('Baseline tri-linear', ...
           'Piecewise-eta model (FOM)', ...
           'POD-ROM piecewise-eta (10 modes)', ...
           'Location','NorthEast');
    grid on;
end
