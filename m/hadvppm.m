function [out] = hadvppm( nn, dt, dx, con, vel, mscl, flxarr, flux1, flux2 )

%      include "camx.prm"
%      real con(nn),vel(nn),flxarr(nn),mscl(nn),saflux(nn)
%     real*8 flux1,flux2

TWO3RDS = 2./3.;

fm = zeros(nn,1);
fp = zeros(nn,1);
cm = zeros(nn,1);
dc = zeros(nn,1);

%Zero order polynomial at the boundary cells

cm(2)  = con(2);
cm(nn) = con(nn-1);

%First order polynomial at the next cells, no monotonicity constraint needed

cm(3) = (con(3) + con(2))/2.;
cm(nn-1) = (con(nn-1) + con(nn-2))/2.;


%Second order polynomial inside the domain

for i = 3 : nn-2
    %Compute average slope in the i'th cell
    dc(i) = 0.5*(con(i+1) - con(i-1));
    
    %Guarantee that CM lies between CON(I) and CON(I+1)
    %monotonicity constraint    
    if ((con(i+1) - con(i))*(con(i) - con(i-1))) <= 0.
        dc(i) = sign(dc(i))*min([abs(dc(i)) 2.*abs(con(i+1) - con(i)) 2.*abs(con(i) - con(i-1))]);
    else
        dc(i) = 0.;
    end
end
      
for i = 3 : nn-3
    cm(i+1) = con(i) + 0.5*(con(i+1) - con(i)) + (dc(i) - dc(i+1))/6.;
end

for i = 2 : nn-1
    cr(i) = cm(i+1);
    cl(i) = cm(i);
end

%Generate piecewise parabolic distributions

for i = 2:nn-1
    if (cr(i) - con(i))*(con(i) - cl(i)) >= 0
        dc(i) = cr(i) - cl(i);
        c6(i) = 6.*(con(i) - 0.5*(cl(i) + cr(i)));
        if dc(i)*c6(i) > dc(i)*dc(i)
            cl(i) = 3.*con(i) - 2.*cr(i);
        elseif -dc(i)*dc(i) > dc(i)*c6(i)
            cr(i) = 3.*con(i) - 2.*cl(i);
        end
    else
        cl(i) = con(i)
        cr(i) = con(i)
    end
    dc(i) = cr(i) - cl(i)
    c6(i) = 6.*(con(i) - 0.5*(cl(i) + cr(i)))
end

% Compute fluxes from the parabolic distribution

for i = 2:nn-1
    x = max(0., -vel(i-1)*(dt/dx));
    fm(i) = x*(cl(i) + 0.5*x*(dc(i) + c6(i)*(1. - TWO3RDS*x)));
    %if (x.ge.1) write(*,*)'Courant number is bigger than 1',x ! jgj 10/6/06
    x = max(0., vel(i)*(dt/dx));
    %if (x.ge.1) write(*,*)'Courant number is bigger than 1',x ! jgj 10/6/06
    fp(i) = x*(cr(i) - 0.5*x*(dc(i) - c6(i)*(1. - TWO3RDS*x)));
end

%Compute fluxes from boundary cells assuming uniform distribution

if vel(1) > 0
    x = vel(1)*(dt/dx)
    fp(1) = x*con(1)
end

if vel(nn-1) < 0
    x = -vel(nn-1)*(dt/dx);
    fm(nn) = x*con(nn);
end

%Update concentrations

flxarr(1) = (fp(1) - fm(2))*(dx/dt)
for i = 2 : nn-1
    flxarr(i) = (fp(i) - fm(i+1))*(dx/dt)
    con(i) = con(i) - mscl(i)*(flxarr(i) - flxarr(i-1))*(dt/dx)
end
flux1 = mscl(2)*flxarr(1)
flux2 = mscl(nn-1)*flxarr(nn-1)
end
