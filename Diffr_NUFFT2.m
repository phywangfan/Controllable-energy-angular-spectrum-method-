%% diffractioin process in 2D using NUFFT
Hn = exp(1i*k*(z*sqrt(1-(fxn.*lam).^2-(fyn.*lam).^2)));
Hn=reshape(Hn,Nfx*Nfy,1);
fxn=reshape(fxn,Nfx*Nfy,1);
fyn=reshape(fyn,Nfx*Nfy,1);
AS_nu = nufft2d3(Nx*Ny,x,y,reshape(sig,Nx*Ny,1),-1,eps,Nfx*Nfy,(fxn)*Lx,(fyn)*Ly); 
FH=AS_nu.*Hn; % diffractioin process
E = nufft2d3(Nfx*Nfy,(fxn)*Lx,(fyn)*Ly,FH,1,eps,Nx*Ny,x,y); % diffractive lighfiled 
E=reshape(E,[Ny,Nx]);
E = E./max(abs(E(:)));
% E = E(Nx/2-Nx/4+1:Nx/2+Nx/4,Ny/2-Ny/4+1:Ny/2+Ny/4);
time=toc;
