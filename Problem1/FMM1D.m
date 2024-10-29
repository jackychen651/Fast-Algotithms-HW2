function u = FMM1D(x, q)
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
% A1_interval(j,l), A2_interval(j,l) denotes multipole expansion in the 
% j-th interval at level l
% B1_interval(j,l), B2_interval(j,l) denotes local expansion in the j-th 
% interval at level l

P = size(x, 1);
L = round(0.5*log2(P)); % TODO: change later
M = 2^L;
x_idx = zeros(P, L);
A1_interval = zeros(2^L, L);
A2_interval = zeros(2^L, L);
B0_interval = zeros(2^L, L);
B1_interval = zeros(2^L, L);
B2_interval = zeros(2^L, L);
mid = zeros(2^L, L);
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
            mid(j,l)=-1+(2*j-1)/2^l;
            A1_interval(j,l) = sum(q(I{j,l}),1)/2;
            A2_interval(j,l) = -dot(x(I{j,l})-mid(j,l), q(I{j,l}))/2;
        end
    else
        for j=1:2^l
            mid(j,l)=-1+(2*j-1)/2^l;
            A1_interval(j,l)=A1_interval(2*j-1,l+1)+A1_interval(2*j,l+1);
            A2_interval(j,l)=A2_interval(2*j-1,l+1)-(mid(2*j-1,l+1)-mid(j,l))*A1_interval(2*j-1,l+1)...
                            +A2_interval(2*j,l+1)-(mid(2*j,l+1)-mid(j,l))*A1_interval(2*j,l+1);
        end
    end
end
for l=2:L
    for i = 1:2^l
        B0_interval(i,l)=B0_interval(ceil(i/2),l-1);
        B1_interval(i,l)=B1_interval(ceil(i/2),l-1);
        B2_interval(i,l)=B2_interval(ceil(i/2),l-1);
        for j=1:2^l
            
            if abs(j-i)<=1 || abs(ceil(j/2) - ceil(i/2)) >1
                continue
            end
            B0_interval(i,l)=B0_interval(i,l)+abs(mid(j,l)-mid(i,l))*A1_interval(j,l)-sign(j-i)*A2_interval(j,l);
            B2_interval(i,l)=B2_interval(i,l)-sign(j-i)*A1_interval(j,l);
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

for i=1:2^L
    for k=1:size(I{i,L},1)
        u(I{i,L}(k))=u(I{i,L}(k))+B0_interval(i,L)+B1_interval(i,L)*abs(x(I{i,L}(k))-mid(i,L))+B2_interval(i,L)*abs(x(I{i,L}(k))-mid(i,L));
    end
end