
function  [delx,dely,delt,X_dual,Y_dual,H_z_new] = KF_2d_Lorentz(k_mac)
             

%[X_dual_Y_main,Y_main_X_dual,E_x_new]
%[X_main_Y_dual,Y_dual_X_main,E_y_new]
%[X_dual,Y_dual,H_z_new]

mu0=1;       % Free space permeability
eps0 = 1;
c = 1;


ND=20; delx=1/ND;     % Avoid dispersion   
Nx= round(1/delx);
Ny = Nx;
delx = 1/Nx;
dely = delx;

Final_T =  0.5;
delt= (0.9/c)*(1/sqrt((1/delx^2)+(1/dely^2))); 
NT = round((Final_T)/delt); 
delt = (Final_T)/NT;

delx = delx/k_mac;
dely = delx;
delt = delt/k_mac;
NT = NT*k_mac;
Nx = Nx*k_mac;
Ny = Ny*k_mac;

%X's and T's
x = 0:delx:Nx*delx;
y = 0:dely:Ny*dely;
x_dual = x(2:end) - delx/2;
y_dual = y(2:end) - dely/2;

[X_dual,Y_dual] = meshgrid(x_dual,y_dual);
[X_dual_Y_main,Y_main_X_dual] = meshgrid(x_dual,y);
[X_main_Y_dual,Y_dual_X_main] = meshgrid(x,y_dual);

%Initial Data
E_x_old = zeros(size(Y_main_X_dual));
E_y_old = zeros(size(Y_dual_X_main));
H_z_old = zeros(size(X_dual));

J_x_old = zeros(size(Y_main_X_dual));
J_y_old = zeros(size(Y_dual_X_main));
P_x_old = zeros(size(Y_main_X_dual));
P_y_old = zeros(size(Y_dual_X_main));

% %PEC Initial Data
% wave_k = sqrt(2*pi^2);
% theta = give_theta(wave_k,nu);
% E_x_old =  -theta.*exp(delt*theta/2).*cos(pi*X_dual_Y_main).*sin(pi*Y_main_X_dual);
% E_y_old =  theta.*exp(delt*theta/2).*sin(pi*X_main_Y_dual).*cos(pi*Y_dual_X_main);
% H_z_old = ((wave_k^2)./pi).*cos(pi*X_dual).*cos(pi*Y_dual);
% 
% J_x_old =  -(theta^2+wave_k^2).*exp(delt*theta/2).*cos(pi*X_dual_Y_main).*sin(pi*Y_main_X_dual);
% J_y_old = (theta^2+wave_k^2).*exp(delt*theta/2).*sin(pi*X_main_Y_dual).*cos(pi*Y_dual_X_main);
% 
% P_x_old = -(-theta - wave_k^2/theta).*exp(delt*theta/2).*cos(pi*X_dual_Y_main).*sin(pi*Y_main_X_dual);
% P_y_old = (-theta - wave_k^2/theta).*exp(delt*theta/2).*sin(pi*X_main_Y_dual).*cos(pi*Y_dual_X_main);


%KF scheme parameters
[Ex_eps_inf,Ex_nu,Ex_omega_0,Ex_omega_p,Ex_sigma] = lorentz_parameters(X_dual_Y_main,Y_main_X_dual,delx,dely);

Ex_KF_J1 = 1/delt + Ex_nu + Ex_omega_0.^2*delt./4;
Ex_KF_J2 = 1/delt - Ex_nu - Ex_omega_0.^2*delt./4;
Ex_KF_J = Ex_KF_J2./Ex_KF_J1;
Ex_KF_E1 = 1 + (Ex_sigma*delt)./(2*eps0.*Ex_eps_inf) + (delt*eps0.*Ex_omega_p.^2)./(4*eps0.*Ex_eps_inf.*Ex_KF_J1);
Ex_KF_E2 = 1 - (Ex_sigma*delt)./(2*eps0.*Ex_eps_inf) - (delt*eps0.*Ex_omega_p.^2)./(4*eps0.*Ex_eps_inf.*Ex_KF_J1);
Ex_KF_E = Ex_KF_E2./Ex_KF_E1;

[Ey_eps_inf,Ey_nu,Ey_omega_0,Ey_omega_p,Ey_sigma] = lorentz_parameters(X_main_Y_dual,Y_dual_X_main,delx,dely);

Ey_KF_J1 = 1/delt + Ey_nu + Ey_omega_0.^2*delt./4;
Ey_KF_J2 = 1/delt - Ey_nu - Ey_omega_0.^2*delt./4;
Ey_KF_J = Ey_KF_J2./Ey_KF_J1;
Ey_KF_E1 = 1 + (Ey_sigma*delt)./(2*eps0.*Ey_eps_inf) + (delt*eps0.*Ey_omega_p.^2)./(4*eps0.*Ey_eps_inf.*Ey_KF_J1);
Ey_KF_E2 = 1 - (Ey_sigma*delt)./(2*eps0.*Ey_eps_inf) - (delt*eps0.*Ey_omega_p.^2)./(4*eps0.*Ey_eps_inf.*Ey_KF_J1);
Ey_KF_E = Ey_KF_E2./Ey_KF_E1;

%Source
source = b((delt/2):delt:(Final_T-delt/2),4);
src_loc = round(size(X_dual,1)*4/10);


%     vidfile = VideoWriter('hmm_vid.avi');
%     vidfile.FrameRate = 60;
%     open(vidfile);

for i = 1:NT
    if(mod(i,round(NT/100)) == 0)
        disp(round(i/round(NT/100)))

        % figure(1)
        % imagesc(x_dual,y_dual,H_z_old,[-1 1])
        % pause(0.01)
    end



    E_x_new = Ex_KF_E.*E_x_old + ...
        (delt./(eps0.*Ex_eps_inf))./Ex_KF_E1.*...
            (matrix_padding(diff(H_z_old,1,1),1)./dely ...
            -(0.5).*(1+Ex_KF_J).*J_x_old + Ex_omega_0.^2./(2*Ex_KF_J1).*P_x_old);
    E_y_new = Ey_KF_E.*E_y_old + ...
        (delt./(eps0.*Ey_eps_inf))./Ey_KF_E1.*...
            (-matrix_padding(diff(H_z_old,1,2),2)./delx ...
             -(0.5).*(1+Ey_KF_J).*J_y_old + Ey_omega_0.^2./(2*Ey_KF_J1).*P_y_old);  

    J_x_new = Ex_KF_J.*J_x_old - Ex_omega_0.^2./Ex_KF_J1.*P_x_old ...
        + eps0.*Ex_omega_p.^2.*(E_x_new + E_x_old)./(2*Ex_KF_J1);
    J_y_new = Ey_KF_J.*J_y_old - Ey_omega_0.^2./Ey_KF_J1.*P_y_old ...
        + eps0.*Ey_omega_p.^2.*(E_y_new + E_y_old)./(2*Ey_KF_J1);

    P_x_new = P_x_old + delt.*(J_x_new + J_x_old)./2;
    P_y_new = P_y_old + delt.*(J_y_new + J_y_old)./2;

    H_z_new = H_z_old - delt./mu0.*(diff(E_y_new,1,2)./delx - diff(E_x_new,1,1)./dely);
    H_z_new(src_loc:src_loc+1,src_loc:src_loc+1) = ...
        H_z_new(src_loc:src_loc+1,src_loc:src_loc+1)...
        + (delt)./(delx*dely)./mu0.*source(i)/2;
    
    %Update variables
    E_x_old = E_x_new;
    E_y_old = E_y_new;
    H_z_old = H_z_new;

    J_x_old = J_x_new;
    J_y_old = J_y_new;

    P_x_old = P_x_new;
    P_y_old = P_y_new;




%         fig = figure('visible','off');
%         plot(x_dual,E_new)
%         xlabel('x')
%         ylabel('Electric Field')
%         title('Numerical Solution with HMM')
%         axis([0 Nx*delx -0.2 1])
%         pause(1/1000)
%         writeVideo(vidfile,getframe(gcf));

end

%          close(vidfile)

end

function out = matrix_padding(in,dim)
    
    if dim == 1
        pad_len = size(in,2);
        out = [zeros(1,pad_len); in; zeros(1,pad_len)];        
    elseif dim == 2
        pad_len = size(in,1);
        out = [zeros(pad_len,1) in zeros(pad_len,1)];   
    end

end

function [eps_inf,nu,omega_0,omega_p,sigma] = lorentz_parameters(X,Y,delx,dely)

    if(any(size(X)~=size(Y)))
        error("Size of X and Y are not the same")
    end

    eps_inf = zeros(size(X));
    nu = zeros(size(X));
    omega_0 = zeros(size(X));
    omega_p = zeros(size(X));
    sigma = zeros(size(X));

    for i = [-1,1]
        for j = [-1,1]
            
            eps_inf(X+i*delx/4<=(0.5)) = eps_inf(X+i*delx/4<=(0.5)) + 2;
            eps_inf(X+i*delx/4>(0.5)) = eps_inf(X+i*delx/4>(0.5)) + 1;

            nu(X+i*delx/4<=(0.5)) = nu(X+i*delx/4<=(0.5)) + 1.25;
            nu(X+i*delx/4>(0.5)) = nu(X+i*delx/4>(0.5)) + 0;

            omega_0(X+i*delx/4<=(0.5)) = omega_0(X+i*delx/4<=(0.5)) + 2*pi;
            omega_0(X+i*delx/4>(0.5)) = omega_0(X+i*delx/4>(0.5)) + 0;

            omega_p(X+i*delx/4<=(0.5)) = omega_p(X+i*delx/4<=(0.5)) + 1;
            omega_p(X+i*delx/4>(0.5)) = omega_p(X+i*delx/4>(0.5)) + 0;

            % sigma(X+i*delx/4<=(0.5)) = sigma(X+i*delx/4<=(0.5)) + 0;
            % sigma(X+i*delx/4>(0.5)) = sigma(X+i*delx/4>(0.5)) +0;

        end
    end

    eps_inf = eps_inf./4;
    nu = nu./4;
    omega_0 = omega_0./4;
    omega_p = omega_p./4;
    sigma = sigma./4;
    
end

function theta = give_theta(k,nu)
    f = @(x) x.^4 - 2*nu*x.^3 + (2+k^2).*x.^2 - 2*nu*k^2.*x + k^2;
    theta = fzero(f,[0,1]);
end


function out = b(t,f_bw)

    an = [0.353222222,-0.488,0.145,-0.010222222];
    out = zeros(size(t));
    ind = t > 0 & t < 1.55/f_bw;

    for i = 1:length(an)
        %out(ind) = out(ind) + an(i).*cos(2*pi*(i-1)*(f_bw/1.55).*t(ind));
        out(ind) = out(ind) - (i-1).*an(i).*sin(2*pi*(i-1)*(f_bw/1.55).*t(ind));
    end

end