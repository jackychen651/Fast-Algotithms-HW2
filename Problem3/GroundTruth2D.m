function u = GroundTruth2D(r_x, theta_x, r_y, theta_y, q_y)

P = size(r_x, 1);
u = zeros(P, 1);
for i = 1:P
    for j = 1:P
        u(i) = u(i) + q_y(j)/(r_x(i)*exp(1j*theta_x(i))-r_y(j)*exp(1j*theta_y(j)));
    end
end