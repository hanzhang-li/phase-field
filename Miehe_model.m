% 本程序用于求解应力应变矩阵C，应力和应变能
% C--Fourth order tensor（四阶张量，应力应变关系）；strain_en--strain energy（应变能）；
% stress_vector 应力向量；  stress_p 拉应力;  phi_p 拉应变能
function [C,stress_p,stress_vector,strain_en,phi_p]=Miehe_model(strain_tensor,phase,lemada,meu,xk)
    V2=[1 1;2 2;1 2];
    heav=@(X)1/2*(X==0)+1*(X>0);
    [eig_vec,eig_val]=eig(strain_tensor);
    eigenvalues_p=eig_val;
    eigenvalues_p(eigenvalues_p<0)=0;
    eigenvalues_n=eig_val-eigenvalues_p;
    trace_strain=eig_val(1,1)+eig_val(2,2);
    J=[1 1 0;1 1 0;0 0 0];
    C=((1-phase)^2+xk)*lemada*heav(trace_strain)*J+lemada*heav(-trace_strain)*J;
    trace_strain_p=trace_strain;
    trace_strain_p(trace_strain_p<0)=0;
    trace_strain_n=trace_strain-trace_strain_p;
    stress_tensor_p=zeros(2,2);
    stress_tensor_n=zeros(2,2);
    phi_p=lemada/2*trace_strain_p^2+meu*(eigenvalues_p(1,1)^2+eigenvalues_p(2,2)^2);
    phi_n=lemada/2*trace_strain_n^2+meu*(eigenvalues_n(1,1)^2+eigenvalues_n(2,2)^2);
    strain_en=((1-phase)^2+xk)*phi_p+phi_n;
    for s=1:2
        stress_tensor_p=stress_tensor_p+2*meu*eigenvalues_p(s,s)*eig_vec(:,s)*eig_vec(:,s)';
        stress_tensor_n=stress_tensor_n+2*meu*eigenvalues_n(s,s)*eig_vec(:,s)*eig_vec(:,s)';
    end
    stress_tensor_p=stress_tensor_p+lemada*trace_strain_p*eye(2);
    stress_tensor_n=stress_tensor_n+lemada*trace_strain_n*eye(2);
    stress_tensor=((1-phase)^2+xk)*stress_tensor_p+stress_tensor_n;
    stress_p=[stress_tensor_p(1,1) stress_tensor_p(2,2) stress_tensor_p(1,2)];
    stress_vector=[stress_tensor(1,1) stress_tensor(2,2) stress_tensor(1,2)];
    for m=1:3
        for n=1:3
            for a=1:2
                HEA_P=heav(eig_val(a,a));
                HEA_N=heav(-eig_val(a,a));
                TP_1=2*meu*HEA_P*eig_vec(V2(m,1),a)*eig_vec(V2(m,2),a)*eig_vec(V2(n,1),a)*eig_vec(V2(n,2),a);
                TN_1=2*meu*HEA_N*eig_vec(V2(m,1),a)*eig_vec(V2(m,2),a)*eig_vec(V2(n,1),a)*eig_vec(V2(n,2),a);
                C(m,n)=C(m,n)+((1-phase)^2+xk)*TP_1+TN_1;
                for b=1:2
                    if(eig_val(a,a)~=eig_val(b,b)&&a~=b)
                        TEN1=(eig_vec(V2(m,1),a)*eig_vec(V2(m,2),b)*eig_vec(V2(n,1),a)*eig_vec(V2(n,2),b));
                        TEN2=(eig_vec(V2(m,1),a)*eig_vec(V2(m,2),b)*eig_vec(V2(n,1),b)*eig_vec(V2(n,2),a));
                        TP_2=meu*(eigenvalues_p(a,a)-eigenvalues_p(b,b))/(eig_val(a,a)-eig_val(b,b))*(TEN1+TEN2);
                        TN_2=meu*(eigenvalues_n(a,a)-eigenvalues_n(b,b))/(eig_val(a,a)-eig_val(b,b))*(TEN1+TEN2);
                        C(m,n)=C(m,n)+((1-phase)^2+xk)*TP_2+TN_2;
                    elseif(eig_val(a,a)==eig_val(b,b)&&a~=b)
                        HEA_P=heav(eig_val(a,a));
                        HEA_N=heav(-eig_val(a,a));
                        TEN1=(eig_vec(V2(m,1),a)*eig_vec(V2(m,2),b)*eig_vec(V2(n,1),a)*eig_vec(V2(n,2),b));
                        TEN2=(eig_vec(V2(m,1),a)*eig_vec(V2(m,2),b)*eig_vec(V2(n,1),b)*eig_vec(V2(n,2),a));
                        TP_3=meu*HEA_P*(TEN1+TEN2);
                        TN_3=meu*HEA_N*(TEN1+TEN2);
                        C(m,n)=C(m,n)+((1-phase)^2+xk)*TP_3+TN_3;                  
                    end
                end
            end
        end
    end
    C=(C+C')/2;
end