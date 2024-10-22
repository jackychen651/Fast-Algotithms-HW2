function u = GroundTruth1D(x, q)
% Notations and Assumptions
% u is interaction, x is point coordinates, q is charges
% Assume we have P charges, and we divide the intervals to M parts
% We also assume that all points are within [-1, 1]

% Therefore, the j-th interval is [-1+2*(j-1)/M,-1+2*j/M].
% We introduce x_idx to be the index of the interval that x is in. For
% example, for x[j], its index will be ceil((x(j)+1)/(2/M)).
% We also introduce I as a cell. I[j] is the  indices of all points
% in [-1+2*(j-1)/M,-1+2*j/M]

P = size(x, 1);
u = zeros(P, 1);
for i = 1:P
    for j = 1:P
        u(i) = u(i) + q(j) * abs(x(i) - x(j)) / 2;
    end
end