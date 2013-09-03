nn = 5;
dt = 900;
dx = 36000;
conInit = zeros(nn,1) + 1;
velInit = zeros(nn,1) + 0.5;
mscl = zeros(nn,1) + 1;
flxarr = zeros(nn,1);

flux1 = 0;
flux2 = 0;

[con, flux1, flux2, fluxarr, fm, fp, cm, cl, cr, dc, c6] = hadvppm(nn, dt, dx, conInit, velInit, mscl, flxarr, flux1, flux2);