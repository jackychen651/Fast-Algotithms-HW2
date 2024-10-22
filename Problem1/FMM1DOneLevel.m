function u = FMM1DOneLevel(x, q)
% Notations and Assumptions
% u is interaction, x is point coordinates, q is charges
% Assume we have P charges, and we divide the intervals to M parts
% We also assume that all points are within [-1, 1]
% The decision of M is unimportant in this problem. We set M to be sqrt(P)
% Therefore, the j-th interval is [-1+2*(j-1)/M,-1+2*j/M].
% We introduce x_idx to be the index of the interval that x is in. For
% example, for x[j], its index will be ceil((x(j)+1)/(2/M)).
% We also introduce I as a cell. I[j] is the  indices of all points
% in [-1+2*(j-1)/M,-1+2*j/M]

P = size(x, 1);
M = round(P^0.5);
x_idx = ceil((x+1)/(2/M));
I = cell(1, M);
u = zeros(P, 1);
for i = 1:P
    I{x_idx(i)} = [I{x_idx(i)}; i];
end
for i=1:M
    for j=1:M
        if abs(j-i)<=1
            continue
        end
        x_c = -1+(2*j-1)/M;
        hat_q_0 = sum(q(I{j}), 1);
        hat_q_1 = dot(x(I{j}), q(I{j})) - x_c * hat_q_0;
        for k=1:size(I{i},1)
            u(I{i}(k)) = u(I{i}(k)) + hat_q_0 * abs(x(I{i}(k))-x_c) / 2 - hat_q_1 * sign(x(I{i}(k)) - x_c) / 2;
        end
    end
    for j=max(1, i - 1):min(M, i + 1)
        for k=1:size(I{i},1)
            for l=1:size(I{j},1)
                u(I{i}(k)) = u(I{i}(k)) + q(I{j}(l))*abs(x(I{i}(k))-x(I{j}(l))) / 2;
            end
        end
    end
end
