nn = 150;
dt = 1.0;
dx = 36000;
x = 1:1:nn;
sigma = 1;
%conInit = zeros(nn,1) + 1;
conInit = exp(-(x-10).^2 ./ (2*sigma^2));
velInit = zeros(nn,1) + 150;
mscl = zeros(nn,1) + 1;
flxarr = zeros(nn,1);

flux1 = 0;
flux2 = 0;
con = conInit;
vel = velInit;
for i = 1:20000
    
    [con flxarr flux1 flux2] = hadvppm(nn, dt, dx, con, vel, mscl, flxarr, flux1, flux2);
    
end

plot(x,conInit,x,con)
legend('Initial Concentration', 'Final Concentration')
title('hadvppm.m Gaussian wave packet test case')