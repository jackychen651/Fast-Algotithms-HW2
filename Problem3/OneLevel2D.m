function u = OneLevel2D(r, theta)

z=r*exp(1j*theta);

P = size(r_x, 1);
M = round(P^0.5);
p=10; %change later
x_idx = ceil((Real(x)+1)/(2/M));
y_idx = ceil((Real(y)+1)/(2/M));
for i_x=1:M
    for i_y=1:M
        if abs(i_x-i_y)<=1 
            pass
        for j_x=1:M
            for j_y=1:M

            end
        end
        end
    