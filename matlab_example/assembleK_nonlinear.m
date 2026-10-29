% =========================================================
%                   서브함수들
% =========================================================

function [K, reg] = assembleK_nonlinear(umax, N, uy1, uy2, k1, k2, k3)
% assembleK_nonlinear :
%   각 스토리마다 "지금까지 경험한 최대 상대변위" umax(i)를 기준으로
%   3구간 선형 강성(k1, k2, k3)을 결정하고, 전체 전단건물 강성행렬 K를 조립.
%   reg(i) = 1 -> k1,  2 -> k2,  3 -> k3
%
%   입력:
%     umax(i) : i층 스토리가 시간이력 중 경험한 상대변위의 최대 절대값
%
%   한 번 uy1, uy2를 넘으면 강성은 회복되지 않음 (누적 손상 모델).

K   = zeros(N);
reg = zeros(N,1);   % 각 스토리의 현재 구간 (1,2,3)

for i = 1:N
    
    du = umax(i);   % 이미 abs()가 적용된 최대 상대변위라고 가정
    
    % 3구간 선형 강성 및 구간 index
    if du <= uy1
        ke      = k1;
        reg(i)  = 1;
    elseif du <= uy2
        ke      = k2;
        reg(i)  = 2;
    else
        ke      = k3;
        reg(i)  = 3;
    end
    
    % 스토리 스프링을 2자유도 트러스처럼 조립
    if i == 1
        K(1,1) = K(1,1) + ke;
    else
        K(i,  i  ) = K(i,  i  ) + ke;
        K(i-1,i-1) = K(i-1,i-1) + ke;
        K(i,  i-1) = K(i,  i-1) - ke;
        K(i-1,i  ) = K(i-1,i  ) - ke;
    end
end

end
