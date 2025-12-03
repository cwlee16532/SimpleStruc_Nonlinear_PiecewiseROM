%% ===== 서브함수: 선형 스토리강성 행렬 조립 =====
function K = assembleK_linear(N, k_story)
% 모든 층에 동일한 스토리 강성 k_story를 가진 전단건물의
% 선형 강성행렬 K를 조립

K = zeros(N);

for i = 1:N
    if i == 1
        K(1,1) = K(1,1) + k_story;
    else
        K(i,  i  ) = K(i,  i  ) + k_story;
        K(i-1,i-1) = K(i-1,i-1) + k_story;
        K(i,  i-1) = K(i,  i-1) - k_story;
        K(i-1,i  ) = K(i-1,i  ) - k_story;
    end
end

end