clear;

% outside constants
lambda_c = 1550e-9;



% simulation required
n_p = 2.4507;
n_m = 2.4296;

n1 = 1.9;
n2 = 1.8;

theta_taper1 = 0;
theta_taper2 = 0;

E0 = 1;
E1 = E0;
E2 = 0;

dn = n_p - n_m;

L_num = 20;
L1 = linspace(4e-6, 20e-6, L_num);
L2 = linspace(1e-6, 6e-6, L_num);


%t = cos(pi*dn/lambda*L1);
%k = sin(pi*dn/lambda*L1);

C = @(l1) [cos(pi*dn/lambda_c*l1) 1i*sin(pi*dn/lambda_c*l1) ; 1i*sin(pi*dn/lambda_c*l1) cos(pi*dn/lambda_c*l1)]*exp(1i*pi/lambda_c*(n_p - n_m));
Pt = [exp(1i*theta_taper1) 0;0 exp(1i*theta_taper2)];
P = @(l2) [exp(1i*2*n1*pi/lambda_c*l2) 0; 0  exp(1i*2*n2*pi/lambda_c*l2)];

Eout = zeros(L_num);
E3 = zeros(L_num);
E4 = zeros(L_num);
eta_bar = zeros(L_num);
eta_cross = zeros(L_num);
for i = 1:L_num
    for j = 1:L_num
        Eout = C(L1(i)) * Pt * P(L2(j)) * Pt * C(L1(i));
        E3(i,j) = Eout(1);
        E4(i,j) = Eout(2);

        I3 = abs(Eout(1))^2;
        I4 = abs(Eout(2))^2;
        Iout = I3 + I4;
        eta_bar(i,j) = I3/Iout;
        eta_cross(i,j) = I4/Iout;
    end
end

figure;
[L1,L2] = meshgrid(L1,L2);
pcolor(L1, L2, eta_bar);
hold on
[C,h] = contour(L1, L2, eta_bar, [0.5 0.5], 'r', 'LineWidth', 3);
colorbar
%hold off






