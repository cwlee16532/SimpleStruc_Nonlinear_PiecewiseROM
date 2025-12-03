function K = assembleK_from_story_stiffness(N, k_story_vec)
% k_story_vec(i) : i번째 스토리 스프링 강성
K = zeros(N);

for i = 1:N
    ke = k_story_vec(i);
    if i == 1
        % ground - 1층
        K(1,1) = K(1,1) + ke;
    else
        K(i,i)       = K(i,i)       + ke;
        K(i-1,i-1)   = K(i-1,i-1)   + ke;
        K(i,i-1)     = K(i,i-1)     - ke;
        K(i-1,i)     = K(i-1,i)     - ke;
    end
end

end
