%Arbitrary announcement date T_ann, applied to:
%Ireland (2007) NK Model, replicates Cagliarini and Kulish (2013, Fig 3)
%To study a different example, simply change the parameters and matrices
%Model structures are defined in the 'Insert' files

clc; clear; close all;

% Announcement date and final date before terminal structure
T_ann = 4; T_tild = 7;
T_sim = 15; %Simulation length

% Model and calibration
%run Insert_NK_forward_guidance
run Insert_NK_inflation_target

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
 
 %Prepare for simulations
 X = X_init;
 
 %For comparison
 X_orig = X; X_fin = X_init_new; 

%Simulation results        
for t=1:T_sim 
        
        X = Omeg_t(:,:,t)*X + Gama_t(:,:,t)*e_vec(:,t) + Psi_t(:,:,t);
        
        %Store for later
        X_stack(:,t) = X;
        
        %Under original structure
        X_orig = Omega_bar*X_orig + Gama_bar*e_vec(:,t) + Psi_bar;
        X_stack_orig(:,t) = X_orig;
        
        %Under terminal structure
        X_fin = Omega_tild*X_fin + Gama_tild*e_vec(:,t) + Psi_tild;
        X_stack_fin(:,t) = X_fin;
    
        Periods(t) = t;
    
end 

NK_inflation_target_plotter
