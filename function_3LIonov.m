function [pt_gen_up,pt_gen_low]=function_3LIonov(wvec,Rhof,Rho_0,Rho_1,Kf,Vp_0,Vp_1,L0,...
                                     B,Porosity,PHI_IONOV,K,...
                                     kpvec_0,kpvec_1,kTvec_0,kTvec_1,Ap0,Ap1);

% Ionov and Maximov (1996) : 10.1111/j.1365-246X.1996.tb05643.x
% Minato et al. (2021) : arXiv:2112.03410 [physics.geo-ph]
% 3-layer model (Top layer and Bottom layer have the same properties)
% Exact solution

varrho__f=Rhof;
varrho__1=Rho_1;
varrho__0=Rho_0;
K__f=Kf;
V__p0=Vp_0;
V__p1=Vp_1;
L=L0;
varphi=Porosity;

Z=0;
I=i;
nw=length(wvec);
pt_gen_up=zeros(1,nw);
pt_gen_low=zeros(1,nw);
for iw=1:nw
    w=wvec(iw);
    k__p0=kpvec_0(iw);
    k__p1=kpvec_1(iw);
    k__0=kTvec_0(iw);
    k__1=kTvec_1(iw);
    A__p0=Ap0(iw);
    A__p1=Ap1(iw);
    D__e0=1/(i*k__p0);
    Phi=PHI_IONOV(iw);

    

    pt_gen_up(iw)=-16*varrho__f*(...
               -(1/4)*(k__p0*varrho__1-k__p1*varrho__0)*(k__1+k__p1)*(-A__p0*V__p1^2*k__p0*varrho__0*K__f*k__p1^2+V__p1^2*(k__p0*A__p0*(varrho__0-varrho__1)*k__1+varrho__0*A__p1*(-k__0+k__p0)*(k__0+k__p0))*K__f*k__p1+A__p0*V__p1^2*k__p0*varrho__1*K__f*k__1^2+B*varphi*varrho__0*K*w^2*Phi*(-k__0+k__p0)*(k__0+k__p0))*(k__0+k__1)*exp(-I*L*(k__1-k__p1)+I*k__0*Z)...
               +(1/4)*(-k__1+k__0)*(k__1+k__p1)*(-A__p0*V__p1^2*k__p0*varrho__0*K__f*k__p1^2+V__p1^2*(k__p0*A__p0*(varrho__0-varrho__1)*k__1+varrho__0*A__p1*(-k__0+k__p0)*(k__0+k__p0))*K__f*k__p1+A__p0*V__p1^2*k__p0*varrho__1*K__f*k__1^2+B*varphi*varrho__0*K*w^2*Phi*(-k__0+k__p0)*(k__0+k__p0))*(k__p0*varrho__1+k__p1*varrho__0)*exp(I*L*(k__1-k__p1)+I*k__0*Z)...
               -(1/4)*(-k__1+k__p1)*(k__0+k__1)*(k__p0*varrho__1+k__p1*varrho__0)*(-A__p0*V__p1^2*k__p0*varrho__0*K__f*k__p1^2-V__p1^2*(k__p0*A__p0*(varrho__0-varrho__1)*k__1-varrho__0*A__p1*(-k__0+k__p0)*(k__0+k__p0))*K__f*k__p1+A__p0*V__p1^2*k__p0*varrho__1*K__f*k__1^2+B*varphi*varrho__0*K*w^2*Phi*(-k__0+k__p0)*(k__0+k__p0))*exp(-I*L*(k__1+k__p1)+I*k__0*Z)...
               +(1/4)*(-k__1+k__0)*(-k__1+k__p1)*(k__p0*varrho__1-k__p1*varrho__0)*(-A__p0*V__p1^2*k__p0*varrho__0*K__f*k__p1^2-V__p1^2*(k__p0*A__p0*(varrho__0-varrho__1)*k__1-varrho__0*A__p1*(-k__0+k__p0)*(k__0+k__p0))*K__f*k__p1+A__p0*V__p1^2*k__p0*varrho__1*K__f*k__1^2+B*varphi*varrho__0*K*w^2*Phi*(-k__0+k__p0)*(k__0+k__p0))*exp(I*L*(k__1+k__p1)+I*k__0*Z)...
               +exp(I*k__0*Z)*varrho__0*(-A__p0*V__p1^2*k__p0*varrho__1*K__f*k__p1^2+V__p1^2*K__f*A__p1*(k__0+k__p0)*(-k__0*varrho__0+k__p0*varrho__1)*k__p1+A__p0*V__p1^2*k__p0*varrho__1*K__f*k__1^2+B*varphi*K*w^2*Phi*(k__0+k__p0)*(-k__0*varrho__0+k__p0*varrho__1))*k__p1*k__1*(-k__0+k__p0))*k__p0*w^2*D__e0...
        /(((k__p0*varrho__1+k__p1*varrho__0)^2*exp(-I*k__p1*L)-exp(I*k__p1*L)*(k__p0*varrho__1-k__p1*varrho__0)^2)*(k__0+k__p0)*(-k__1+k__p1)*V__p1^2*(k__1+k__p1)*(-k__0+k__p0)*K__f*((k__0+k__1)^2*exp(-I*k__1*L)-exp(I*k__1*L)*(-k__1+k__0)^2));


    pt_gen_low(iw)=...
        -8*varrho__f*(-(k__p0*varrho__1+k__p1*varrho__0)*k__1*((-A__p0*V__p1^2*k__p0*K__f*k__p1^3+A__p1*V__p1^2*K__f*(-k__0+k__p0)*(k__0+k__p0)*k__p1^2+(A__p0*V__p1^2*K__f*k__1^2*k__p0+(-k__0+k__p0)*(k__0+k__p0)*(B*K*Phi*varphi*w^2+A__p1*K__f*V__p1^2*k__0))*k__p1+B*varphi*K*w^2*Phi*k__0*(-k__0+k__p0)*(k__0+k__p0))*varrho__0-A__p0*V__p1^2*k__p0*varrho__1*K__f*k__0*(-k__1+k__p1)*(k__1+k__p1))*exp(I*k__0*Z-I*k__p1*L)...
              +(k__1*(-k__0+k__p0)*(k__0+k__p0)*(B*K*Phi*varphi*w^2+A__p1*K__f*V__p1^2*k__p1)*varrho__0+varrho__1*(-V__p1^2*K__f*A__p0*(k__1+k__p0)*k__p1^2+V__p1^2*K__f*A__p1*(-k__0+k__p0)*(k__0+k__p0)*k__p1+A__p0*K__f*V__p1^2*k__1^3+A__p0*V__p1^2*K__f*k__1^2*k__p0+B*varphi*K*w^2*Phi*(-k__0+k__p0)*(k__0+k__p0))*k__p0)*varrho__0*k__p1*(k__0+k__1)*exp(I*k__0*Z-I*k__1*L)...
              -((-A__p0*V__p1^2*k__p0*K__f*k__p1^3+A__p1*V__p1^2*K__f*(-k__0+k__p0)*(k__0+k__p0)*k__p1^2+(A__p0*V__p1^2*K__f*k__1^2*k__p0+(-k__0+k__p0)*(k__0+k__p0)*(B*K*Phi*varphi*w^2-A__p1*K__f*V__p1^2*k__0))*k__p1-B*varphi*K*w^2*Phi*k__0*(-k__0+k__p0)*(k__0+k__p0))*varrho__0+A__p0*V__p1^2*k__p0*varrho__1*K__f*k__0*(-k__1+k__p1)*(k__1+k__p1))*(k__p0*varrho__1-k__p1*varrho__0)*k__1*exp(I*k__0*Z+I*k__p1*L)...
              -varrho__0*(-k__1+k__0)*k__p1*((-k__1*(-k__0+k__p0)*(k__0+k__p0)*(B*K*Phi*varphi*w^2+A__p1*K__f*V__p1^2*k__p1)*varrho__0+varrho__1*(-V__p1^2*K__f*A__p0*(-k__1+k__p0)*k__p1^2+V__p1^2*K__f*A__p1*(-k__0+k__p0)*(k__0+k__p0)*k__p1-A__p0*K__f*V__p1^2*k__1^3+A__p0*V__p1^2*K__f*k__1^2*k__p0+B*varphi*K*w^2*Phi*(-k__0+k__p0)*(k__0+k__p0))*k__p0)*exp(I*k__0*Z+I*k__1*L)...
                                             ))*k__p0*w^2*D__e0...
    /(((k__0+k__1)^2*exp(-I*k__1*L)-exp(I*k__1*L)*(-k__1+k__0)^2)*(k__0+k__p0)*(k__1+k__p1)*V__p1^2*(-k__1+k__p1)*((k__p0*varrho__1+k__p1*varrho__0)^2*exp(-I*k__p1*L)-exp(I*k__p1*L)*(k__p0*varrho__1-k__p1*varrho__0)^2)*(-k__0+k__p0)*K__f);

end




