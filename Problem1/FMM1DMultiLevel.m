function u = FMM1DMultiLevel(x, q)
% Notations and Assumptions
% u is interaction, x is point coordinates, q is charges
% L is total level, l is current level
% Assume we have P charges, and we divide the intervals to M parts
% M=2^L
% We also assume that all points are within [-1, 1]
% The j-th interval in level l is is [-1+2*(j-1)/2^l,-1+2*j/2^l].
% We introduce x_idx to be a cell, containg all indices of the 
% (level, interval) pair that x is in. For example, for x[j], its index 
% will be (l, ceil((x(j)+1)/(2/2^l))). To simplify, we store it as column 
% vector. i.e., x_idx(j, l) = ceil((x(j)+1)/(2/2^l)), 1<=j<=P, 1<=l<=L
% We also introduce I as a cell. I[j,l] is the indices of all points in the
% j-th interval at level l. i.e., [-1+2*(j-1)/2^l,-1+2*j/2^l].
% 1<=j<=P, 1<=l<=L
% q0_interval(j,l) denotes the sum of charges in the j-th interval at level 
% l
% q1_interval(j,l) denotes the sum of weighted charges in the j-th interval
% at level l. i.e., q1_interval(j,l) = \sum_k q_k*x_k

P = size(x, 1);
L = round(0.65*log2(P)); % TODO: change later
M = 2^L;
x_idx = zeros(P, L);
q0_interval = zeros(2^L, L);
for j=1:P
    for l=1:L
        x_idx(j,l)=ceil((x(j)+1)/(2/2^l));
    end
end
I = cell(M, L);
u = zeros(P, 1);
for l = 1:L
    for j = 1:P
        I{x_idx(j,l),l} = [I{x_idx(j,l),l}; j];
    end
end
for l = L:-1:1
    if l == L
        for j = 1:2^l
            q0_interval(j,l) = sum(q(I{j,l}),1);
            q1_interval(j,l) = dot(x(I{j,l}), q(I{j,l}));
        end
    else
        for j = 1:2^l
            q0_interval(j,l) = q0_interval(2*j-1, l+1) + q0_interval(2*j, l+1);
            q1_interval(j,l) = q1_interval(2*j-1, l+1) + q1_interval(2*j, l+1);
        end
    end
end
for l=2:L
    for i = 1:2^l
        for j=1:2^l
            if abs(j-i)<=1 || abs(ceil(j/2) - ceil(i/2)) >1
                continue
            end
            x_c = -1+(2*j-1)/2^l;
            hat_q_0 = q0_interval(j,l);
            hat_q_1 = q1_interval(j,l) - x_c * hat_q_0;
            for k=1:size(I{i,l},1)
                u(I{i,l}(k)) = u(I{i,l}(k)) + hat_q_0 * abs(x(I{i,l}(k))-x_c) / 2 - hat_q_1 * sign(x(I{i,l}(k)) - x_c) / 2;
            end
        end
    end
end
for i=1:2^L
    for j=max(1, i - 1):min(2^L, i + 1)
        for k=1:size(I{i,L},1)
            for m=1:size(I{j,L},1)
                u(I{i,L}(k)) = u(I{i,L}(k)) + q(I{j, L}(m))*abs(x(I{i, L}(k))-x(I{j,L}(m))) / 2;
            end
        end
    end
end

