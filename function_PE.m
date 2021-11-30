function [PE]=function_PE(PE)
%Calculating poroelastic parameters (including the Skempton
%coefficient) using given vp_dry,vs_dry,rhos,phi, Ks and Kf

PE.rho_dry=(1-PE.phi)*PE.rhos;
PE.rho_bulk=(1-PE.phi)*PE.rhos+PE.phi*PE.rhof;
PE.G=PE.vs_dry^2*PE.rho_dry;
PE.K_dry=PE.vp_dry^2*PE.rho_dry-4/3*PE.G;
PE.alpha=1-PE.K_dry./PE.Ks;
PE.M=( (PE.alpha-PE.phi)./PE.Ks+PE.phi./PE.Kf ).^(-1);
PE.C=PE.M.*PE.alpha;
PE.H=PE.K_dry+4/3*PE.G+PE.M.*PE.alpha.^2; %H=K_undrained+4/3*G (see 2.6c, Guan2011)


%Gassmann K_undrained=alpha^2*M+K_dry (see Guan2011)
PE.Ksat=PE.K_dry+(1-PE.K_dry/PE.Ks)^2/(...
    PE.phi/PE.Kf+(1-PE.phi)/PE.Ks-PE.K_dry/PE.Ks^2);
PE.K_undrained=PE.K_dry+PE.alpha^2*PE.M;
PE.Vp0=sqrt(PE.H/PE.rho_bulk); %==sqrt((K_undrained+4/3G)/rho_bulk)
PE.Vs0=sqrt(PE.G/PE.rho_bulk); %==sqrt((K_undrained+4/3G)/rho_bulk)
%Skempton Coefficient
PE.B=(1/PE.K_dry-1/PE.Ks)/(1/PE.K_dry-1/PE.Ks+PE.phi*(1/PE.Kf-1/PE.Ks));
