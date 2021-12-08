%Modeling_VSP_tubewave_Porous_Layer
% S. Minato
% Delft University of Technology
% s.minato-1@tudelft.nl
%
% This code calculates pressure waveforms at a borehole due to
% a normally incident plane P wave (vertical seismic profiling, VSP).
% The borehole is a circular column filled with fluid, 
% and it is embedded in horizontally layered media. 
%


%Note about the Skempton coefficient:
%      : Flag_Skempton switches the following two conditions.
%      : Flag_Skempton=0 -> pext=-1/3(txx+tyy+tzz).
%      : Flag_Skempton=1 -> pext=-1/3(txx+tyy+tzz)*B.
%      : "B" is the Skempton coefficient. 
%      : "B" is either calculated from given Vp and Vs (not recommended) or PE.B

%===========================================================================
% This version uses modeling parameters designed for a simulation of
% a thin porous layer sandwiched between two elastic half-spaces and
% a borehole with a constant radius
%===========================================================================


clear all
close all

Flag_Skempton=1; 


if(Flag_Skempton)
    disp('==========================')
    disp('pext=-1/3tii*B activated!')
    disp('==========================')
end


%===Construction of a layered model========
%   N-layer model (1 and N are half-spaces)
%==========================================
n=5; % Total number of layers
Vpvec=[5000 5000 5000 5000 5000]; % size n. P-wave velocity at each layer (m/s)
Vsvec=[3000 3000 3000 3000 3000]; % size n. S-wave velocity at each layer (m/s)
Rhovec=[2500 2500 2500 2500 2500]; %Bulk density at each layer (kg/m3)
zn_org=[10 19.5 20.5 30]; % size n-1. The depth of boundaries (m). 
rn=0.055*ones(1,n); %size n. Borehole radius (m)


% Fluid properties
nu_dyn=1E-3; %Fluid dynamic viscosity (Pa*s)
Vf=1500;
Rhof=1000;
Kf=Rhof*Vf^2;%Fluid Bulk modulus

% Porous-layer properties 
Phivec=0.3*ones(1,n); %size n. Porosity (m3/m3)
Kappa0_vec=zeros(1,n); %size n. Static permeability (m^2), 1 Darcy=9.869E-13.
Kappa0_vec(3)=1*9.869E-13; %Here assuming that the 3rd layer is a porous layer

% Skempton coefficient
PE.vp_dry=5170; %dry vp
PE.vs_dry=3198; %dry vs
PE.rhos=3143; %grain density
PE.rhof=Rhof; %fluid density
PE.phi=0.3; %porosity
PE.Ks=10E10; %grain bulk modulus
PE.Kf=Kf;
PE=function_PE(PE);

%Because the 3rd layer is a porour layer, replace the vp and vs at
%this layer with those of the low-frequency limits (Gassmann)
disp('==================================')
disp('Vp,Vs,Rho at porous media replaced!')
disp('==================================')
Vpvec(3)=PE.Vp0;
Vsvec(3)=PE.Vs0;
Rhovec(3)=PE.rho_bulk;



%--------------------------------------
% Waveform modeling parameters
%--------------------------------------
zvec_rec=[0:0.2:45]; %Receiver depth locations (m)

dt=0.25E-4; %sampling interval (sec)
ns=16001; %sample number

tvec=[0:dt:(ns-1)*dt];%time vector
dw=1/tvec(end)*2*pi; %frequency sampling interval (radian)
wvec_org=[0:dw:(ns-1)*dw];

nw_proc=300; %sample number in frequency (< ns/2)
wvec_proc=wvec_org(1:nw_proc); %frequency vector (radian) where
                               %wavefield is calculated

%Note: Here I define the Ricker wavelet for the source signature. 
%    : The propagator matrix method for an elastic wavefield considers 
%    : this source signature as tzz. See below in the 1D propagator matrix method.

% Creating a source wavelet
f0=200; % center frequency (Hz)
delay=1/f0*2; % delay time (s)


Ricker=2*pi^2*f0^2*(1.0 - 2.0*(pi^2)*(f0^2)*((tvec-delay).^2)) .* exp(-(pi^2)*(f0^2)*((tvec-delay).^2)); %Ricker wavelet

Ricker=Ricker./max(abs(Ricker));

fRicker=conj(fft(Ricker)); %MATLAB FT -> Aki-Richards FT
figure;plot(abs(fft(Ricker)))


%Plotting the model
figure;
for in=1:n
    if(in==1)
        subplot(161),hold on;plot(Vpvec(in)*ones(1,2),[zn_org(in)-10 zn_org(in)],'k:','Linewidth',2)
        subplot(162),hold on;plot(Vsvec(in)*ones(1,2),[zn_org(in)-10 zn_org(in)],'k:','Linewidth',2)
        subplot(163),hold on;plot(Rhovec(in)*ones(1,2)*1E-3,[zn_org(in)-10 zn_org(in)],'k:','Linewidth',2)
        subplot(164),hold on;plot(rn(in)*ones(1,2)*1E+3,[zn_org(in)-10 zn_org(in)],'k:','Linewidth',2)
        subplot(165),hold on;plot(Phivec(in)*ones(1,2),[zn_org(in)-10 zn_org(in)],'k:','Linewidth',2)
        subplot(166),hold on;plot(Kappa0_vec(in)*ones(1,2)/9.869E-13,[zn_org(in)-10 zn_org(in)],'k:','Linewidth',2)
    elseif(in==n)
        subplot(161),hold on;plot(Vpvec(in)*ones(1,2),[zn_org(in-1) zn_org(in-1)+10],'k:','Linewidth',2)
        grid on;set(gca,'Ydir','Reverse');title('Vp(m/s)');ylabel('Depth(m)')
        subplot(162),hold on;plot(Vsvec(in)*ones(1,2),[zn_org(in-1) zn_org(in-1)+10],'k:','Linewidth',2)    
        grid on;set(gca,'Ydir','Reverse');title('Vs(m/s)');ylabel('Depth(m)')
        subplot(163),hold on;plot(Rhovec(in)*ones(1,2)*1E-3,[zn_org(in-1) zn_org(in-1)+10],'k:','Linewidth',2)    
        grid on;set(gca,'Ydir','Reverse');title('Rho(g/cc)');ylabel('Depth(m)')
        subplot(164),hold on;plot(rn(in)*ones(1,2)*1E3,[zn_org(in-1) zn_org(in-1)+10],'k:','Linewidth',2)    
        grid on;set(gca,'Ydir','Reverse');title('Radius(mm)');ylabel('Depth(m)')
        subplot(165),hold on;plot(Phivec(in)*ones(1,2),[zn_org(in-1) zn_org(in-1)+10],'k:','Linewidth',2)    
        grid on;set(gca,'Ydir','Reverse');title('\phi');ylabel('Depth(m)')
        subplot(166),hold on;plot(Kappa0_vec(in)*ones(1,2)/9.869E-13,[zn_org(in-1) zn_org(in-1)+10],'k:','Linewidth',2)    
        grid on;set(gca,'Ydir','Reverse');title('\kappa_0 (D)');ylabel('Depth(m)')
    else
        subplot(161),hold on;plot(Vpvec(in)*ones(1,2),[zn_org(in-1) zn_org(in)],'k-','Linewidth',2)
        subplot(161),hold on;plot([Vpvec(in-1) Vpvec(in)],ones(1,2)*zn_org(in-1),'k.-','Linewidth',2)
        subplot(162),hold on;plot(Vsvec(in)*ones(1,2),[zn_org(in-1) zn_org(in)],'k-','Linewidth',2)
        subplot(162),hold on;plot([Vsvec(in-1) Vsvec(in)],ones(1,2)*zn_org(in-1),'k.-','Linewidth',2)
        subplot(163),hold on;plot(Rhovec(in)*ones(1,2)*1E-3,[zn_org(in-1) zn_org(in)],'k-','Linewidth',2)
        subplot(163),hold on;plot([Rhovec(in-1) Rhovec(in)]*1E-3,ones(1,2)*zn_org(in-1),'k.-','Linewidth',2)
        subplot(164),hold on;plot(rn(in)*ones(1,2)*1E3,[zn_org(in-1) zn_org(in)],'k-','Linewidth',2)
        subplot(164),hold on;plot([rn(in-1) rn(in)]*1E3,ones(1,2)*zn_org(in-1),'k.-','Linewidth',2)
        subplot(165),hold on;plot(Phivec(in)*ones(1,2),[zn_org(in-1) zn_org(in)],'k-','Linewidth',2)
        subplot(165),hold on;plot([Phivec(in-1) Phivec(in)],ones(1,2)*zn_org(in-1),'k.-','Linewidth',2)
        subplot(166),hold on;plot(Kappa0_vec(in)*ones(1,2)/9.869E-13,[zn_org(in-1) zn_org(in)],'k-','Linewidth',2)
        subplot(166),hold on;plot([Kappa0_vec(in-1) Kappa0_vec(in)]/9.869E-13,ones(1,2)*zn_org(in-1),'k.-','Linewidth',2)
    end
    
end


%%===SHIFTING===
shift_z=0;
zn_proc=zn_org+shift_z; %
%%===SHIFTING===

%================================================================
%   1D Propagator-matrix method: Calculation of elastic wavefield
%   due to a normally incident plane P wave
%================================================================

uEvec_allfreq=zeros(2,n,nw_proc);%Upgoing and Downgoing potential amplitudes [Un;Dn] for each layer and all frequencies. "E" stands for Elastic wavefield.

for iw=2:nw_proc
    
    w=wvec_proc(iw);
    kp=w./Vpvec;
    
    %Constructing M matrix. size (2,2) and 1 to n-1
    Mn=zeros(2,2,n-1); 
    for in=1:n-1
        kp1=kp(in);
        kp2=kp(in+1);
        z1=zn_proc(in);
        rho1=Rhovec(in);
        rho2=Rhovec(in+1);
        
        a1=rho2/rho1;
        a2=kp2/kp1;
        
        m11=1.0/(2.0)*(a1+a2)*exp(i*(kp1-kp2)*z1); %note: Aki-Richards FT
        m12=1.0/(2.0)*(a1-a2)*exp(i*(kp1+kp2)*z1); 
        m21=1.0/(2.0)*(a1-a2)*exp(-i*(kp1+kp2)*z1); 
        m22=1.0/(2.0)*(a1+a2)*exp(-i*(kp1-kp2)*z1); 
        
        Mn(:,:,in)=[m11 m12; m21 m22];
    end


    %Calculating MT matrix, (2,2)
    MT=eye(2,2);
    for in=1:n-1
        MT=MT*Mn(:,:,in);    
    end
    %------------------------------
    % Solution at each layer
    %------------------------------
    uEvec=zeros(2,n);%=[Un;Dn] and each layer


    % Radiation condition at the top and the bottom layers
    
    % Definition of 'unit' amplitude source can be defined here.
    % Unit displacement -> (1/(i*kp(1))
    % Unit partivle velocity -> (1/(iw*i*kp(1))) 
    % Unit stress -> (-Rhovec(1)*w^2)
    D1=1/(-Rhovec(1)*w^2); %this simulate unit amplitudes in tzz, and defined source wavelet is applied in this unit
    
    %Additional phase shift due to shift_z
    D1=D1*exp(i*kp(1)*(-shift_z)); 

    Dn=D1/MT(2,2);
    uEvec(:,n)=[0.0;Dn];

    % Solving the potentials by connecting from below
    for in=n-1:-1:1
        uEvec(:,in)=Mn(:,:,in)*uEvec(:,in+1);    
    end

    % Check if uvec(:,1)==[U1;0]
    U1=MT(1,2)/MT(2,2)*D1;
    if (abs(uEvec(1,1)-U1) > 1E-5*abs(U1))
        disp("Warning:U1 is not retrieved after connecting Haskel matrix from below??")
        tmp_str=sprintf('iw=%d,U1=%f,uEvec(1,1)=%f\n',iw,U1,uEvec(1,1));
        disp(tmp_str)
    end

    uEvec_allfreq(:,:,iw)=uEvec;

end %freq


% $$$ return

%%Calculating time-domain waveforms using the potentials (Elastic wave)

fdata=zeros(length(zvec_rec),ns);

for irec=1:length(zvec_rec)
    znow=zvec_rec(irec)+shift_z;
    
    %Detecting a layer number
    tmp=zn_proc-znow;
    tmp(find(tmp<0))=Inf;
    [tmp1 tmp2]=min(tmp(:));
    in_now=tmp2;
    if(sum(isinf(tmp))==length(tmp))
        in_now=length(zn_proc)+1; %last layer
    end
    
    U_amp=squeeze(uEvec_allfreq(1,in_now,:)); %upgoing,layer in_now
    D_amp=squeeze(uEvec_allfreq(2,in_now,:)); %downgoing,layer in_now

    kpvec=wvec_proc/Vpvec(in_now);    
    
    phi=U_amp.'.*exp(-i*kpvec*znow)+D_amp.'.*exp(i*kpvec*znow); %displacement potential (note: kc homogenous)
    
    
    tmpdata=-i*kpvec.*U_amp.'.*exp(-i*kpvec*znow)...
            +i*kpvec.*D_amp.'.*exp(i*kpvec*znow); %displacement (uz)
    tmpdata=-i*wvec_proc.*tmpdata; %particle velocity (vz)
    

    tmpdata=tmpdata.*fRicker(1:nw_proc); %Note: Unit incident wave is
                                         % assumed in potential
                                         % (see the phase difference between phi and tmpdata above)
    fdata(irec,1:nw_proc)=tmpdata;
end

fdata(:,1)=0;
fdata=conj(fdata); %Aki-Richards FT -> MATLAB FT
data=ifft(fdata,[],2,'symmetric');
data_E=data.';

figure;imagesc(zvec_rec,tvec,data_E);colorbar;
title("Elastic wave (VZ)")
xlabel('Receiver depth (m)');ylabel('Time (s)')


%==================================================================
%   1D Propagator-matrix method: Calculation of a borehole pressure 
%   field using the elastic-wavefield potential amplitudes
%================================================================

%Necessary elastic moduli and velocities
Evec=Rhovec.*Vsvec.^2.*(3*Vpvec.^2-4*Vsvec.^2)./(Vpvec.^2-Vsvec.^2); %size n. Young's modulus
mu_vec=Rhovec.*Vsvec.^2; %size n. Shear modulus
lambda_vec=Rhovec.*Vpvec.^2-2*mu_vec; %size n. Lame constant
Kvec=lambda_vec+2/3*mu_vec; %size n. Bulk modulus
CT0vec=sqrt(Vf^2./(1+Rhof*Vf^2./(Rhovec.*Vsvec.^2))); %size n. Tube-wave velocity without porous-layer effects

%Porous-layer effects: Diffusivity (Ionov, eq. A3)
Diffvec=sqrt(Kappa0_vec*Kf./(nu_dyn*Phivec));
TFvec=rn.^2./Diffvec.^2; %(See Ionov after eq. A6) Note: kappa0=0 becomes TF=Inf. In this case, there are no porous-layer effects.



%
uvec_allfreq=zeros(2,n,nw_proc);%Upgoing and Downgoing tube-wave potential amplitudes [Un;Dn] for each layer and all frequencies.

Cy2vec=Evec/Rhof;


if(Flag_Skempton)
    disp('==========================')
    disp('pext=-1/3tii*B activated!')
    disp('==========================')
end

for iw=2:nw_proc
    
    w=wvec_proc(iw);
    kp=w./Vpvec;
    
    %Porous layer effects: Calculating the function PHI of Ionov
    %   and correponding tube-wave velocities
    %Note: all variables related to the porous layer assumes MATLAB-FT
    tmp=sqrt(i*w*TFvec); %be carefull INF when kappa0=0
    PHI_IONOV=tmp.^(-1).*besselk(1,tmp)./besselk(0,tmp); %be careful NaN when kappa0=0

    %Tube-wave velocity including the porous-layer effects
    CTvec=CT0vec.*sqrt(...
          1./(1+2*Phivec.*Rhof/Kf.*CT0vec.^2.*PHI_IONOV));

    %Removing NaN when kappa0=0 (those CTs are CT0)
    CTvec(find(Kappa0_vec==0))=CT0vec(find(Kappa0_vec==0));
    CTvec=conj(CTvec); %MATLAB FT -> Aki-Richards FT
    
    kn=w./CTvec;   

    %Constructing M matrix. size (2,2) and 1 to n-1
    Mn=zeros(2,2,n-1); 
    for in=1:n-1
        r1=rn(in);
        r2=rn(in+1);
        k1=kn(in);
        k2=kn(in+1);
        z1=zn_proc(in);
 
        a1=r1^2*k1+r2^2*k2;
        a2=r1^2*k1-r2^2*k2;
        
        m11=a1*exp(i*(k1-k2)*z1); %note: Aki-Richards FT
        m12=a2*exp(i*(k1+k2)*z1); 
        m21=a2*exp(-i*(k1+k2)*z1); 
        m22=a1*exp(-i*(k1-k2)*z1); 

        Mn(:,:,in)=1/(2*r1^2*k1)*[m11 m12; m21 m22];
    end

    %Constructing S matrix. size (2,1) and 1 to n-1
    Sn=zeros(2,n-1); 
    for in=1:n-1        
        r1=rn(in);
        r2=rn(in+1);
        k1=kn(in);
        z1=zn_proc(in);
        
        kp1=kp(in);
        Cy_square=Cy2vec(in);
        E1=Evec(in);
        Vp_now=Vpvec(in);
        Vs_now=Vsvec(in);
        CT_now=CTvec(in);
        Porosity_now=Phivec(in);
        PHI_IONOV_now=conj(PHI_IONOV(in)); %MATLAB FT -> Aki-Richards FT
        K_now=Kvec(in); %Bulk modulus
        
        %Evalauting the discontinuities (delta_p and delta_vz)
        %due to external effective stress.
        %I am at z1=zn(in). Thus, z0=zn(in-1)
        if (in==1) %In this case, z0 is -Inf
           I1=0;
           I2=0;
           I3=0;
           I4=0;
        else
            z0=zn_proc(in-1);

            %Integartion element I1:
            I1=exp(i*k1*z1)/i/(-k1+kp1)*(exp(i*(-k1+kp1)*z1)-exp(i*(-k1+kp1)*z0))...
               -exp(-i*k1*z1)/i/(k1+kp1)*(exp(i*(k1+kp1)*z1)-exp(i*(k1+kp1)*z0));
            %Integartion element I2:
            I2=exp(i*k1*z1)/i/(-k1-kp1)*(exp(i*(-k1-kp1)*z1)-exp(i*(-k1-kp1)*z0))...
               -exp(-i*k1*z1)/i/(k1-kp1)*(exp(i*(k1-kp1)*z1)-exp(i*(k1-kp1)*z0));
            %Integartion element I3:
            I3=exp(i*k1*z1)/i/(-k1+kp1)*(exp(i*(-k1+kp1)*z1)-exp(i*(-k1+kp1)*z0))...
               +exp(-i*k1*z1)/i/(k1+kp1)*(exp(i*(k1+kp1)*z1)-exp(i*(k1+kp1)*z0));
            %Integartion element I4:
            I4=exp(i*k1*z1)/i/(-k1-kp1)*(exp(i*(-k1-kp1)*z1)-exp(i*(-k1-kp1)*z0))...
               +exp(-i*k1*z1)/i/(k1-kp1)*(exp(i*(k1-kp1)*z1)-exp(i*(k1-kp1)*z0));

        end        
        %---the discontinuities (delta_p and delta_vz) due to elastic wave---
        U_Eamp=squeeze(uEvec_allfreq(1,in,iw)); %upgoing elastic wave,layer in
        D_Eamp=squeeze(uEvec_allfreq(2,in,iw)); %downgoing elastic wave,layer in

        A1P=w^2/kp1*(1/(2*Vs_now^2)-1/Vp_now^2);
        delta_p_A=-i*w*Rhof*CT_now*kp1*A1P...
                  *(D_Eamp*I1+U_Eamp*I2);
        delta_vz_A=-i*w*kp1*A1P...
                  *(D_Eamp*I3+U_Eamp*I4);
        
        
 

        if (Kappa0_vec(in)~=0)
            %the discontinuities (delta_p and delta_vz) due to a porous layer 
            %(PHI_IONOV_now from Aki-Richards FT)
            delta_p_B=Rhof*CT_now*(-i*w*Porosity_now/Kf*PHI_IONOV_now)*w^2/Vp_now^2*K_now...
                      *(D_Eamp*I1+U_Eamp*I2);
            delta_vz_B=(-i*w*Porosity_now/Kf*PHI_IONOV_now)*w^2/Vp_now^2*K_now...
                *(D_Eamp*I3+U_Eamp*I4);
            
        
            %%Option B
            if(Flag_Skempton)
                % The following two lines calculates the Skempton coefficient (B)
                % assuming the given Vp and Vs to be those in the drained condition.
                % This is not recommended.
% $$$                 K_d=K_now; %assuming Vp and Vs are drained
% $$$                 B=1-1/(1+Kf/(Porosity_now*K_d)); %skepmton coefficient
                % The following line uses the Skempton coefficient (B) 
                % provided in PE.B (recommended)
                B=PE.B;

                %Application of the Skempton coefficient (B)
                delta_p_B=B*delta_p_B;
                delta_vz_B=B*delta_vz_B;                
            end
            
        else
            delta_p_B=0;
            delta_vz_B=0;
        end

 
        %Calculating Sn=-An^(-1)*[delta_p;rn^2*delta_vz]
        Aninv=[exp(i*k1*z1)/(2*w^2*Rhof) -exp(i*k1*z1)/(2*w*pi*r1^2*k1);
               exp(-i*k1*z1)/(2*w^2*Rhof) exp(-i*k1*z1)/(2*w*pi*r1^2*k1)]; 

    
        %Additional injection source (q) due to the step-like change in the radius
        tmpdata=-i*kp1.*U_Eamp.'.*exp(-i*kp1*z1)...
                +i*kp1.*D_Eamp.'.*exp(i*kp1*z1); %displacement (uz)
        tmpdata=-i*w.*tmpdata; %particle velocity (vz)
        
        %---Kurkjian's dV due to P-wave vz:velocity---
        dVn=pi*(r2^2-r1^2)*tmpdata;  %dV/dt=pi(r2^2-r1^2)*vz
        dvz=dVn/(pi*r1^2); %displacement src = qsrc
                           %=dV/dt/(pi*r1^2)@media 1
        %-----------------------------------

        %The total discontinuities
        dp_total=delta_p_A+delta_p_B;
        dvz_total=delta_vz_A+delta_vz_B+dvz;
        
        %The source term
        Sn(:,in)=1/(2*Rhof*w^2*k1)*[(Rhof*w*dvz_total-k1*dp_total)*exp(i*k1*z1);...
                                   -(Rhof*w*dvz_total+k1*dp_total)*exp(-i*k1*z1)];


    end

    
    %Calculating MT matrix, (2,2), and ST matrix, (2,1)
    MT=eye(2,2);
    ST=zeros(2,1);
    for in=1:n-1
        ST=ST+MT*Sn(:,in);
        MT=MT*Mn(:,:,in);    
    end

    %------------------------------
    % Solution at each layer
    %------------------------------
    uvec=zeros(2,n);%=[Un;Dn] and each layer

     
    % Radiation condition at the top and the bottom layers
            
    % Bottom layer Un is known from
    % the amplitude of incidence wave at the bottom layer. 
    % This assumes that there is a downgoing wave only in the elastic wavefield.
    D_Eamp=squeeze(uEvec_allfreq(2,n,iw)); %downgoing elastic wave,layer in
    
    z0=zn_proc(n-1);
    k1=kn(n);
    kp1=kp(n);
    Vp_now=Vpvec(n);
    Vs_now=Vsvec(n);
    CT_now=CTvec(n);

    %Integration element I1:
    I1=2*k1/i/(kp1^2-k1^2)*(exp(i*kp1*z0));

    %Integration element I3:
    I3=2*kp1/i/(kp1^2-k1^2)*(exp(i*kp1*z0));


    A1P=w^2/kp1*(1/(2*Vs_now^2)-1/Vp_now^2);
    delta_p_A=-i*w*Rhof*CT_now*kp1*A1P...
              *(D_Eamp*I1);
    delta_vz_A=-i*w*kp1*A1P...
        *(D_Eamp*I3);
    
    Un=exp(i*k1*z0)/(2*Rhof*w^2*k1)*(...
        k1*delta_p_A-Rhof*w*delta_vz_A...
        );
    
    

    % Top layer D1 is known from
    % the amplitude of incidence wave at the top layer.
    % This assumes downgoing and upgoing waves in the elastic wavefield.
    D_Eamp=squeeze(uEvec_allfreq(2,1,iw)); %downgoing elastic wave,layer in
    U_Eamp=squeeze(uEvec_allfreq(1,1,iw)); %downgoing elastic wave,layer in
    
    %z0 is -inf
    z1=zn_proc(1);
    k1=kn(1);
    kp1=kp(1);
    Vp_now=Vpvec(1);
    Vs_now=Vsvec(1);
    CT_now=CTvec(1);

    %Integration element I1:
    I1=2*k1/i/(kp1^2-k1^2)*(exp(i*kp1*z1));

    %Integration element I2:
    I2=2*k1/i/(kp1^2-k1^2)*(exp(-i*kp1*z1));
    
    %Integration element I3:
    I3=2*kp1/i/(kp1^2-k1^2)*(exp(i*kp1*z1));
    
    %Integration element I4:
    I4=-2*kp1/i/(kp1^2-k1^2)*(exp(-i*kp1*z1));

    A1P=w^2/kp1*(1/(2*Vs_now^2)-1/Vp_now^2);
    delta_p_A=-i*w*Rhof*CT_now*kp1*A1P...
              *(D_Eamp*I1+U_Eamp*I2);
    delta_vz_A=-i*w*kp1*A1P...
        *(D_Eamp*I3+U_Eamp*I4);
    
    
    D1=exp(-i*k1*z1)/(2*Rhof*w^2*k1)*(...
        k1*delta_p_A+Rhof*w*delta_vz_A...
        );
    
    % Setting up the relation [U1;D1]=MT*[Un;Dn]+ST
    Dn=(D1-MT(2,1)*Un-ST(2))/MT(2,2);
    U1=MT(1,1)*Un+MT(1,2)/MT(2,2)*(D1-MT(2,1)*Un-ST(2))+ST(1);
    
    uvec(:,n)=[Un;Dn];

        
    % Solving the potentials by connecting from below
    for in=n-1:-1:1
        uvec(:,in)=Mn(:,:,in)*uvec(:,in+1)+Sn(:,in);    
    end

    uvec_allfreq(:,:,iw)=uvec;

end %freq


%%Calculating time-domain waveforms using the potentials (Borehole)

fdata=zeros(length(zvec_rec),ns);

for irec=1:length(zvec_rec)
    znow=zvec_rec(irec)+shift_z;
    
    %Detecting a layer number
    tmp=zn_proc-znow;
    tmp(find(tmp<0))=Inf;
    [tmp1 tmp2]=min(tmp(:));
    in_now=tmp2;
    if(sum(isinf(tmp))==length(tmp))
        in_now=length(zn_proc)+1; %last layer
    end

    
    U_amp=squeeze(uvec_allfreq(1,in_now,:)); %upgoing,layer in_now
    D_amp=squeeze(uvec_allfreq(2,in_now,:)); %downgoing,layer in_now


    if(Kappa0_vec(in_now)~=0)
        
        %Porous layer effects: Calculating the function PHI of Ionov
        %   and correponding tube-wave velocities
        %Note: all variables related to the porous layer assumes MATLAB-FT
        tmp=sqrt(i*wvec_proc*TFvec(in_now)); %be carefull INF when kappa0=0
        PHI_IONOV=tmp.^(-1).*besselk(1,tmp)./besselk(0,tmp); %be careful NaN when kappa0=0

        %Tube-wave velocity including the porous-layer effects
        CTvec=CT0vec(in_now).*sqrt(...
            1./(1+2*Phivec(in_now)*Rhof/Kf.*CT0vec(in_now).^2.*PHI_IONOV)); %this is a frequency-vector version.

        
        CTvec=conj(CTvec); %MATLAB FT -> Aki-Richards FT
        
        kcvec=wvec_proc./CTvec;
    else     
        kcvec=wvec_proc./CT0vec(in_now);
    end
    

    phi=U_amp.'.*exp(-i*kcvec*znow)+D_amp.'.*exp(i*kcvec*znow); %displacement potential (note: kc homogenous)
    tmpdata=Rhof*wvec_proc.^2.*phi; %pressure
                                 
    
    %Evalauting the discontinuities (delta_p and delta_vz)
    %due to external effective stress.
    %Note:multiple-freq vector version    
    %I am at z1=zvec(irec). z0 is the nearest boundary in the negative z direction
    
    
    if(in_now==1)
        z1=znow; 
        z0=zn_proc(1);
    elseif(in_now==n)
        z1=znow; 
        z0=zn_proc(n-1);        
    else
        z1=znow; 
        z0=zn_proc(in_now-1);
    end    
    
    %I am at in=in_now
    k1vec=kcvec;    
    kp1vec=wvec_proc/Vpvec(in_now);
    Cy_square=Cy2vec(in_now);
    E1=Evec(in_now);
    Vp_now=Vpvec(in_now);
    Vs_now=Vsvec(in_now);


    %Integration element I1:
    I1=exp(i*k1vec.*z1)/i./(-k1vec+kp1vec).*(exp(i*(-k1vec+kp1vec)*z1)-exp(i*(-k1vec+kp1vec)*z0))...
       -exp(-i*k1vec.*z1)/i./(k1vec+kp1vec).*(exp(i*(k1vec+kp1vec)*z1)-exp(i*(k1vec+kp1vec)*z0));
    %Integration element I2:
    I2=exp(i*k1vec.*z1)/i./(-k1vec-kp1vec).*(exp(i*(-k1vec-kp1vec)*z1)-exp(i*(-k1vec-kp1vec)*z0))...
       -exp(-i*k1vec*z1)/i./(k1vec-kp1vec).*(exp(i*(k1vec-kp1vec)*z1)-exp(i*(k1vec-kp1vec)*z0));
    %Integration element I3:
    I3=exp(i*k1vec*z1)/i./(-k1vec+kp1vec).*(exp(i*(-k1vec+kp1vec)*z1)-exp(i*(-k1vec+kp1vec)*z0))...
       +exp(-i*k1vec*z1)/i./(k1vec+kp1vec).*(exp(i*(k1vec+kp1vec)*z1)-exp(i*(k1vec+kp1vec)*z0));
    %Integration element I4:
    I4=exp(i*k1vec*z1)/i./(-k1vec-kp1vec).*(exp(i*(-k1vec-kp1vec)*z1)-exp(i*(-k1vec-kp1vec)*z0))...
       +exp(-i*k1vec*z1)/i./(k1vec-kp1vec).*(exp(i*(k1vec-kp1vec)*z1)-exp(i*(k1vec-kp1vec)*z0));
    
    if(Kappa0_vec(in_now)~=0)
        CT_now=CTvec; %vector
    else
        CT_now=CT0vec(in_now); %scalar
    end            

    %---the discontinuities (delta_p and delta_vz) due to elastic wave---
    U_Eamp=squeeze(uEvec_allfreq(1,in_now,:)); %upgoing elastic wave,layer in_now
    D_Eamp=squeeze(uEvec_allfreq(2,in_now,:)); %downgoing elastic wave,layer in_now
    
    U_Eamp=U_Eamp.';
    D_Eamp=D_Eamp.';
    
    A1P=wvec_proc.^2./kp1vec*(1/(2*Vs_now^2)-1/Vp_now^2);
    delta_p_A=-i*wvec_proc*Rhof.*CT_now.*kp1vec.*A1P...
              .*(D_Eamp.*I1+U_Eamp.*I2);

    %the discontinuities (delta_p and delta_vz) due to a porous layer 
    %(PHI_IONOV_now from Aki-Richards FT)
    if(Kappa0_vec(in_now)~=0)        
        Porosity_now=Phivec(in_now);                
        K_now=Kvec(in_now);

        delta_p_B=Rhof*CT_now.*(-i*wvec_proc*Porosity_now/Kf.*conj(PHI_IONOV)).*wvec_proc.^2/Vp_now^2*K_now...
                  .*(D_Eamp.*I1+U_Eamp.*I2); %BUGFIX
    
        %%Option B
        if(Flag_Skempton)
                % The following two lines calculates the Skempton coefficient (B)
                % assuming the given Vp and Vs to be those in the drained condition.
                % This is not recommended.
% $$$                 K_d=K_now; %assuming Vp and Vs are drained
% $$$                 B=1-1/(1+Kf/(Porosity_now*K_d)); %skepmton coefficient
                % The following line uses the Skempton coefficient (B) 
                % provided in PE.B (recommended)
                B=PE.B;

                %Application of the Skempton coefficient (B)
                delta_p_B=B*delta_p_B;
        end
    
    else     
        delta_p_B=0;
    end
    
    
    
    tmpdata=tmpdata+delta_p_A+delta_p_B;
    tmpdata=tmpdata.*fRicker(1:nw_proc); %Note: "Unit incident wave" is
                                         % defined in elastic-wave potentials

    
    fdata(irec,1:nw_proc)=tmpdata;
end

    
fdata(:,1)=0;
fdata=conj(fdata); %Aki-Richards FT -> MATLAB FT
data=ifft(fdata,[],2,'symmetric');
data_B=data.';

figure;imagesc(zvec_rec,tvec,data_B);colorbar;
title("Borehole response")
xlabel('Receiver depth (m)');ylabel('Time (s)')


%return

%%Comparison with analytical solutions
disp('Comparison with analytical solutions...')


in_1=3; %the layer number of a thin porous layer
L0=zn_org(3)-zn_org(2); %porous layer thickness 
in_0=1; %the layer number of an elastic layer 

tmp=sqrt(i*wvec_proc*TFvec(in_1)); %be carefull INF when kappa0=0
PHI_IONOV=tmp.^(-1).*besselk(1,tmp)./besselk(0,tmp); %be careful NaN when kappa0=0
CTvec=CT0vec(in_1).*sqrt(...
    1./(1+2*Phivec(in_1)*Rhof/Kf.*CT0vec(in_1).^2.*PHI_IONOV));

CTvec=conj(CTvec); %MATLAB FT -> Aki-Richards FT. Note: This CTvec
                   %is within a porous layer

CTvec_homo=CT0vec(in_0)*ones(size(wvec_proc)); %Note: This CTvec is within a surrounding elastic media


%Top nd Bottom layers have the same properties
Vp0=Vpvec(in_0);
Vs0=Vsvec(in_0);
Rho0=Rhovec(in_0);

%Middle
Vp1=Vpvec(in_1);
Vs1=Vsvec(in_1);
Rho1=Rhovec(in_1);

%---start
CT0=CTvec_homo;
CT1=CTvec;

kp0=wvec_proc./Vp0;
kp1=wvec_proc./Vp1;
k0=wvec_proc./CT0;
k1=wvec_proc./CT1;

AP0=wvec_proc.^2./kp0*(1/(2*Vs0^2)-1/Vp0^2);
AP1=wvec_proc.^2./kp1*(1/(2*Vs1^2)-1/Vp1^2);

Porosity_now=Phivec(in_1);                
K_now=Kvec(in_1); %bulk modulus


%%Exact solution
[pt_gen_up,pt_gen_low]=function_3LIonov(wvec_proc,Rhof,Rho0,Rho1,Kf,Vp0,Vp1,L0,...
                            PE.B,Porosity_now,conj(PHI_IONOV),K_now,...
                            kp0,kp1,k0,k1,AP0,AP1);

De0=1./(i*kp0); %if uz=1
pt_direct0=-Rhof*CT0.*wvec_proc.*2.*k0.*kp0./(kp0.^2-k0.^2).*AP0.*De0;





%------
test_trace_E=squeeze(data_E(:,1)); %elastic VZ. For incident wave
test_trace=squeeze(data_B(:,1)); %for up
TT=abs(zvec_rec(1))/Vpvec(1)+abs(zvec_rec(1)-zn_org(in_1-1))/CT0vec(1)+0.005; %for up

test_trace2=squeeze(data_B(:,end)); %for down
TT2=abs(zvec_rec(1)-zn_org(in_1-1))/Vpvec(end)+abs(zvec_rec(end)-zn_org(in_1-1))/CT0vec(end)+0.005; %for down



tvec_intp=[0:dt:3/f0];
pinc_est=interp1(tvec,test_trace,tvec_intp);
vzinc_est=interp1(tvec,test_trace_E,tvec_intp);
pT_est_up=interp1(tvec-TT,test_trace,tvec_intp);
pT_est_down=interp1(tvec-TT2,test_trace2,tvec_intp);


figure;
subplot(211),plot(tvec,test_trace)
subplot(212),plot(tvec_intp,pinc_est,'k-')
hold on;plot(tvec_intp,pT_est_up,'r-')
hold on;plot(tvec_intp,pT_est_down,'b-')
legend('inc','up','down')
grid on;


% $$$ return

fvec_est=[0:1:(length(tvec_intp)-1)]/tvec_intp(end);

fpinc_est=fft(pinc_est);
fvzinc_est=fft(vzinc_est);
fpT_est=fft(pT_est_up);
fpT_est2=fft(pT_est_down);

ifmax=6;

figure;
subplot(121),plot(fvec_est(2:ifmax),abs(fpT_est(2:ifmax)./fpinc_est(2:ifmax)),'.-')
hold on;plot(wvec_proc/2/pi,abs(pt_gen_up./pt_direct0),'r-') %exact
title('Pup-Pinc ratio')
legend('Haskel','Exact')
grid on;

subplot(122),plot(fvec_est(2:ifmax),abs(fpT_est2(2:ifmax)./fpinc_est(2:ifmax)),'.-')
hold on;plot(wvec_proc/2/pi,abs(pt_gen_low./pt_direct0),'r-') %exact
legend('Haskel','Exact')
title('Pdown-Pinc ratio')
grid on;




