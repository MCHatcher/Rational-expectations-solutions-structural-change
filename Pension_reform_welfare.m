% Arbitrary announcement date T_ann, applied to Application 3 - Pension reform
% Diamond (1965) model with CES utility; see Fedotenkov (2016,EL) and
% Hatcher(2019, EL) for log utility version of the model

clc; clear;

% Announcement date and final date before terminal structure
T_ann = 4; T_tild = 4;
T_sim = 5000; %Simulation length

Nloop = 20;
stack = linspace(0.001,0.999,Nloop);

for z=1:Nloop

SP_discount = stack(z); %social discount factor

% Model and calibration
run Insert_pension_reform

% Fixed structure solutions (Cho and Moreno 2011, JEDC)
run Cho_and_Moreno

%Indicator variable    
ind = ones(T_sim,1); ind(T_tild+1:T_sim,1) = 0;  %permanent structural change
%Initial values for recursion
Omeg = Omega_tild; Gama = Gama_tild; Psi = Psi_tild;  

%Computation of matrix recursion
 for t=T_sim:-1:1       
            
    B1t = ind(t,1)*B1 + (1-ind(t,1))*B1_tild;
    B2t = ind(t,1)*B2 + (1-ind(t,1))*B2_tild;
    B3t = ind(t,1)*B3 + (1-ind(t,1))*B3_tild;
    B4t = ind(t,1)*B4 + (1-ind(t,1))*B4_tild;
    B5t = ind(t,1)*B5 + (1-ind(t,1))*B5_tild;
         
   Omeg = (B1t - B2t*Omeg) \ B3t; 
   Gama = (B1t - B2t*Omeg) \ B4t; 
   Psi = (B1t - B2t*Omeg) \ (B2t*Psi + B5t);    
        
    if t >= T_tild+1
        Omeg = Omega_tild;
        Gama = Gama_tild;
        Psi = Psi_tild;   
    end
    
    if t <= T_ann-1
        Omeg = Omega_bar;
        Gama = Gama_bar;
        Psi = Psi_bar;
    end
                      
    Omeg_t(:,:,t) = Omeg;
    Gama_t(:,:,t) = Gama;
    Psi_t(:,:,t) = Psi;
    
 end

%Initial values
X_init = [0; 0; 0; 0; 0]; 

%Prepare for simulations
X = X_init;
 
%For comparison
X_u = X_init;  %Unannounced reform

%Simulation results        
for t=1:T_sim 
    
        if t < T_ann
            X = Omega_bar*X + Psi_bar;
        end
        
        if t >= T_ann
            X = Omeg_t(:,:,t)*X + Psi_t(:,:,t);
        end
        
        %Store for later
        X_stack(:,t) = X;
        
        if t <= T_ann
            X_u = Omega_bar*X_u + Psi_bar;
        end
        
        if t > T_ann
            X_u = Omega_tild*X_u + Psi_tild;
        end
        
        Xu_stack(:,t) = X_u;
        
        cy(t) = X_stack(1,t); r(t) = X_stack(2,t); w(t) = X_stack(3,t);
        k(t) = X_stack(4,t); co(t) = X_stack(5,t);
        klev(t) = exp(k(t) + log(kstar));
        
        cyu(t) = Xu_stack(1,t); ru(t) = Xu_stack(2,t); wu(t) = Xu_stack(3,t);
        ku(t) = Xu_stack(4,t); cou(t) = Xu_stack(5,t);
    
end 


U_init = (1/(1-sigma))*(exp(cystar))^(1-sigma) + betta*(exp(costar))^(1-sigma) + (exp(cystar))^(1-sigma)*(cy(T_ann-1) +(1-sigma)/2*cy(T_ann-1)^2) + (exp(costar))^(1-sigma)*(co(T_ann) +(1-sigma)/2*co(T_ann)^2);

for t=T_ann:T_sim-1
    
     U(t) = (1/(1-sigma))*( exp(cystar1)^(1-sigma) + betta*exp(costar1)^(1-sigma) ) + (exp(cystar1))^(1-sigma)*(cy(t) +(1-sigma)/2*cy(t)^2) + betta*(exp(costar1))^(1-sigma)*(co(t+1) +(1-sigma)/2*co(t+1)^2);
     Uu(t) = (1/(1-sigma))*( exp(cystar1)^(1-sigma) + betta*exp(costar1)^(1-sigma) ) + (exp(cystar1))^(1-sigma)*(cyu(t) +(1-sigma)/2*cyu(t)^2) + betta*(exp(costar1))^(1-sigma)*(cou(t+1) +(1-sigma)/2*cou(t+1)^2); 
     
   if t <= T_tild
       U(t) = (1/(1-sigma))*( exp(cystar)^(1-sigma) + betta*exp(costar)^(1-sigma) ) + (exp(cystar))^(1-sigma)*(cy(t) +(1-sigma)/2*cy(t)^2) + betta*(exp(costar))^(1-sigma)*(co(t+1) +(1-sigma)/2*co(t+1)^2);
       Uu(t) = (1/(1-sigma))*( exp(cystar)^(1-sigma) + betta*exp(costar)^(1-sigma) ) + (exp(cystar))^(1-sigma)*(cyu(t) +(1-sigma)/2*cyu(t)^2) + betta*(exp(costar))^(1-sigma)*(cou(t+1) +(1-sigma)/2*cou(t+1)^2);
   end
       
   U1(t) = SP_discount^(t-T_ann)*U(t);
   U1u(t) = SP_discount^(t-T_ann)*Uu(t);
   
end

Welfare(z) = sum(U1) + SP_discount^(-1)*U_init;
Welfareu(z) = sum(U1u) + SP_discount^(-1)*U_init;


    lambda(z) = 100*( (Welfare(z)/Welfareu(z))^(1/(1-sigma)) - 1);
    lambda1(z) = 100*(Welfare(z)-Welfareu(z))/Welfareu(z);
    discount(z) = SP_discount;
    zero(z) = 0;
    
    end


figure(1)
hold on, plot(discount,lambda1), hold on, xlabel('Social discount factor (\gamma)'), ylabel('% change in social welfare'), plot(discount,zero,'--k')

figure(2)
hold on, plot(discount, Welfare-Welfareu), hold on, plot(discount,zero,'--k'), xlabel('Social discount factor (\gamma)'), ylabel('Welfare gain of announced reform')





