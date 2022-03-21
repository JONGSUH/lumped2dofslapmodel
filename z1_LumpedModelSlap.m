function [displ,velo] =  z1_LumpedModelSlap(M,C,K,Initial_gap,v_n0,d_n0,F,tt,AnalType)
%%
global dt
epsilon0 = 5e9; 
depsilon = epsilon0*1; % for v0=4.5
epsilon = epsilon0; % for v0=4.5


displ = zeros(size(M,1),length(tt));
velo = zeros(size(M,1),length(tt));
displ(:,1) = d_n0; velo(:,1) = v_n0;

old_E = (1/2)*transpose(velo(:,1))*M*velo(:,1) + (1/2)*transpose(displ(:,1))*K*displ(:,1);
%%
nu = -1; N = [nu; -nu];
%%
flag=0;
Rhamda_n_k = zeros(size(M,1)-1,1);
tol_convergence = 1e-14; 
tol_gaprate=2.85e-26; 


g_n0 = -diff(d_n0)-Initial_gap;
tlist = 0;

d_n1_k_i = d_n0; 
old_delta_d_n1_k_i = 0;

Index_slave_penetration = find(g_n0>0);
[Jacobian_c,Fc,Tn,Gap_Rate] = z2_Jacobian(d_n1_k_i,d_n0,g_n0,Index_slave_penetration,N,Rhamda_n_k,epsilon);

max_iter=100; max_iter_conv = 3;
iter=0;
iter_conv=0;
ii=1;   
tic
while ~flag
    iter = iter+1;
    R_l = ( (2/dt^2)*M + (1/dt)*C + K/2 )*d_n1_k_i - Fc;
    R_r = (F(:,ii) + F(:,ii+1))/2 + ( (2/dt^2)*M + (1/dt)*C - K/2 )*d_n0 + 2/dt*M*v_n0;

    R = R_l - R_r;
    Jacobian = (2/dt^2)*M + K/2 + Jacobian_c;            
    delta_d_n1_k_i = -Jacobian\R;   
    d_n1_k_i_1 = d_n1_k_i + delta_d_n1_k_i;           
    v_n1_k_i_1 = 2*(d_n1_k_i_1 - d_n0)/dt - v_n0;
    new_E = (1/2)*transpose(v_n1_k_i_1)*M*v_n1_k_i_1 + (1/2)*transpose(d_n1_k_i_1)*K*d_n1_k_i_1;

    new_delta_d_n1_k_i = delta_d_n1_k_i;
    
    flag_convergence = norm(new_delta_d_n1_k_i - old_delta_d_n1_k_i);
    flag_energy = 100*(new_E - old_E)/old_E;

    if abs(flag_convergence)<=tol_convergence || iter>=max_iter
        if isempty(Index_slave_penetration) || (max(Gap_Rate)<=1*tol_gaprate) 

            d_n1 = d_n1_k_i_1; v_n1 = v_n1_k_i_1;
            old_E = new_E;

            displ(:,ii+1) = d_n1; velo(:,ii+1) = v_n1;
            d_n0 = d_n1; v_n0 = v_n1;        
            g_n0 = -diff(d_n0)-Initial_gap;       

            Index_slave_penetration = find(g_n0>0);

            d_n1_k_i = d_n0; 
            [Jacobian_c,Fc,Tn0,Gap_Rate] = z2_Jacobian(d_n1_k_i,d_n0,g_n0,Index_slave_penetration,N,Rhamda_n_k,epsilon);
            
            iter=0;
            iter_conv=0;
            tlist(ii+1) = tlist(ii)+dt;
            ii=ii+1;

            if ii==length(tt)
                break
            end 
        else    
            if iter_conv<=max_iter_conv
                iter=0;
                iter_conv = iter_conv+1;
                d_n1_k_i = d_n1_k_i_1;
                if strcmp(AnalType,'AugLag')
                    Rhamda_n_k(Index_slave_penetration) = Tn;
                else
                    Rhamda_n_k(Index_slave_penetration) = 0;
                end

                [Jacobian_c,Fc,Tn,Gap_Rate] = z2_Jacobian(d_n1_k_i,d_n0,g_n0,Index_slave_penetration,N,Rhamda_n_k,epsilon);
            else
                conv_energy = flag_energy;
                iter=0;
                iter_conv=0;
                d_n1_k_i = d_n0;
                if strcmp(AnalType,'AugLag')
                    Rhamda_n_k(Index_slave_penetration) = Tn0;
                else
                    Rhamda_n_k(Index_slave_penetration) = 0;
                end
                epsilon = epsilon + depsilon*1.0;
                [Jacobian_c,Fc,Tn,Gap_Rate] = z2_Jacobian(d_n1_k_i,d_n0,g_n0,Index_slave_penetration,N,Rhamda_n_k,epsilon);
            end
                
        end
    else
        d_n1_k_i = d_n1_k_i_1;
        [Jacobian_c,Fc,Tn,Gap_Rate] = z2_Jacobian(d_n1_k_i,d_n0,g_n0,Index_slave_penetration,N,Rhamda_n_k,epsilon);
    end

end
  
toc
 %%   Plot results for 2DOF
figure(); plot(tlist(1,1:end),displ(1,1:length(tlist))-Initial_gap/2,'r-*'); hold on;
plot(tlist(1,1:end),displ(2,1:length(tlist))+Initial_gap/2,'b-o')
title('Positions')
hold off

figure(); plot(tlist(1,1:end),velo(1,1:length(tlist)),'r-*'); hold on;
plot(tlist(1,1:end),velo(2,1:length(tlist)),'b-o')
title('Velocities')
hold off

total_energy = zeros(1,length(tlist));
mm =diag(M); m1=mm(1); m2=mm(2);
for ii = 1:size(velo,2)    
    kinetic_e = velo(:,ii)'*M*(velo(:,ii));   
    potential_e = displ(:,ii)'*K*(displ(:,ii));
    total_energy(ii) = kinetic_e + potential_e;
end

figure(); plot(tlist(1,1:end),total_energy); grid on; title('Total Energy')
% % % % figure(); plot(displ(1,1:length(tlist)),displ(2,1:length(tlist)),'r-*');
% % % % title('Mode shapes')
% % % % 
% % % % gt = -Initial_gap + displ(1,1:length(tlist)) - displ(2,1:length(tlist));
% % % % gvt = velo(1,1:length(tlist)) - velo(2,1:length(tlist));
% % % % figure()
% % % % plot(gt,gvt,'k:'); title('Phase Plot')
