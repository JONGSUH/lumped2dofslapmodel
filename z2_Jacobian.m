function [Jacobian_c,Fc,Tn,Gap_Rate] = z2_Jacobian(d_n1_k_i,d_n0,g_n0,Index_slave_penetration,N,Rhamda_n_k,epsilon)
global M dt 

Ic=zeros(length(Index_slave_penetration)*2,1);
Xc=zeros(length(Index_slave_penetration)*2,1);

Ij=zeros(length(Index_slave_penetration)*4,1);
Jj=zeros(length(Index_slave_penetration)*4,1);        
Xj=zeros(length(Index_slave_penetration)*4,1);
ntriplets_c=0;
ntriplets_j=0;
Tn = zeros(length(Index_slave_penetration),1);
Gap_Rate = zeros(length(Index_slave_penetration),1);
if isempty(Index_slave_penetration)
    Jacobian_c = zeros(size(M,1),size(M,2));
    Fc=zeros(size(M,1),1);
    Tn = 0;
    Gap_Rate=0;
else

    for jj=1:length(Index_slave_penetration)
        slave_node = Index_slave_penetration(jj);         
        edof = [slave_node; slave_node+1];
        tem_gap_rate = (d_n1_k_i(slave_node) - d_n1_k_i(slave_node+1) - d_n0(slave_node) + d_n0(slave_node+1))/dt;        
        tem_tn = Rhamda_n_k(Index_slave_penetration(jj)) + epsilon * tem_gap_rate;
        
        Tn(jj)=tem_tn;
        Gap_Rate(jj) = tem_gap_rate;

        tem_g_n0 = g_n0(Index_slave_penetration(jj));
        tn = heaviside(tem_g_n0)*heaviside(tem_tn)*abs(tem_tn);

        fc = tn*N;    
        jacobian_c = ( dirac(tem_g_n0)*heaviside(tem_tn)*abs(tem_tn) + heaviside(tem_g_n0)*heaviside(tem_tn)*epsilon/dt ) * N*transpose(N);

        for krow=1:2
            ntriplets_c = ntriplets_c+1;
            Ic(ntriplets_c) = edof(krow);               
            Xc(ntriplets_c) = fc(krow);            
            for kcol=1:2
                ntriplets_j = ntriplets_j+1;
                Ij(ntriplets_j) = edof(krow);
                Jj(ntriplets_j) = edof(kcol);
                Xj(ntriplets_j) = jacobian_c(krow,kcol);                                
            end                
        end
    end
    Jacobian_c = sparse(Ij,Jj,Xj,size(M,1),size(M,2));
    Fc=sparse(Ic,1,Xc,size(M,1),1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
