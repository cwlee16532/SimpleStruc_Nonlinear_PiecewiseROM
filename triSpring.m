function ke = triSpring(u, uy1, uy2, k1, k2, k3)
% triSpring : 3구간 선형 스프링의 탄성계수 계산
%
%   ke = triSpring(u, uy1, uy2, k1, k2, k3)
%
% 입력:
%   u   : 요소 변위 (상대변위)
%   uy1 : 1차 경계 변위 (|u| <= uy1 -> k1)
%   uy2 : 2차 경계 변위 (uy1 < |u| <= uy2 -> k2, 그 이상 -> k3)
%   k1  : 첫 번째 구간 강성
%   k2  : 두 번째 구간 강성
%   k3  : 세 번째 구간 강성
%
% 출력:
%   ke  : 현재 변위에 대한 탄성계수(강성)
%
% 잔류변위는 0으로 가정하고, loading / unloading 에 따른
% 히스테리시스 없이, 단순히 변위 크기에 따른 구간선형 모델임.

du = abs(u);

if du <= uy1
    ke = k1;
elseif du <= uy2
    ke = k2;
else
    ke = k3;
end

end
