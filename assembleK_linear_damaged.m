function K = assembleK_linear_damaged(N, k_story, damageMask, eta)
% assembleK_linear_damaged :
%   모든 층 기본 스토리강성 k_story 이고,
%   damageMask(i) = 1 인 층은 강성을 eta * k_story 로 줄여서
%   선형 전단건물 강성행렬 K를 조립.

K = zeros(N);

for i = 1:N
    if damageMask(i) == 1
        ke = eta * k_story;   % 손상 층
    else
        ke = k_story;         % 비손상 층
    end
    
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
