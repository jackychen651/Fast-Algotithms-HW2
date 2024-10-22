function u = FMM2D(r_x, theta_x, r_y, theta_y, q_y, M)

% We assume that all x will be in (-3,-1)*(-1,1) and y will be in (1,3)*(-1,1)
% Also q_x is the charges of x, q_y is the charges of y.
N=size(r_x,1);
hat_q=zeros(M+1, 1);
u = zeros(N, 1);
for k = 1:M+1
    for j=1:N
        hat_q(k)=hat_q(k)+q_y(j)*r_y(j)^(k-1)*exp(1j*(k-1)*theta_y(j));
    end
end
for i=1:N
    for k=0:M
        u(i)=u(i)+hat_q(k+1) /(r_x(i)^(k+1)*exp(1j*(k+1)*theta_x(i)));
    end
end
