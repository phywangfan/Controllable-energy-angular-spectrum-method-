% clc  
clear
close all
%% fundamental parameters
eta=0.98;   % the weight of user-defined
lam = 532e-6;	% wavelength
k = 2*pi/lam;	% wave number
Nx = 1024;                
Ny = 1024;	% number of the original signal
dp = 0.001;	% sampling interval in spatial domain (square pixel)
Lx = (Nx)*dp;                  
Ly = (Ny)*dp; % source window size
x = linspace(-Lx/2,Lx/2,Nx)';        
y = linspace(Ly/2,-Ly/2,Ny)'; 
[x,y]=meshgrid(x,y);    % cordinate in spatial domain
zc=2*max(Nx,Ny)*dp^2/lam;   % critical distance for transfer function
zn=15;  % zn times the critical distance
z =zc*zn ;  % propagation distance
% ========================= signal setting =======================
% ----------------------- triangle aperture ---------------------
V=[-0.05+0.1,  0.05+0.1, z;
    0+0.1, -0.05+0.1, z;
    0.1+0.1,  0+0.1, z]; % vertex coordinates of triangle
sig=plot_tri(V',Lx,Ly,Nx,Ny); % triangle as a 2D image 
% ----------------------- rectangle aperture ---------------------
%     rx = Ny/2;          
%     ry=  Ny/2;   % numbers of aperture
%     sig = zeros(Ny,Nx);
% %     sig (Ny/2-ry/2+1:Ny/2+ry/2,Nx/2-rx/2+1:Nx/2+rx/2) = 1;  % source window
% sig(x.^2+y.^2<0.1)=1;
% ----------------------- 2D image ---------------------
%     sig=im2double(imread('USAF951.png'));  % import an image
%     sig=imresize(sig(:,:,1),[Ny,Nx]);
    [ry,rx]=size(sig);  % aperture is the window of image
% -------------------------------------------------------
figure,imshow(sig,[]); title('source signal)')
wx=(rx-1)*dp;
wy=(ry-1)*dp; % aperture size
% ================== propagation parameters ================
Lx = (Nx)*dp;                  
Ly = (Ny)*dp; % source window size

dfx=1/dp/(Nx);
dfy=1/dp/(Ny);    % sampling interval in spatial domain
fx = linspace(-1/2/dp,1/2/dp,Nx)';   
fy = linspace(1/2/dp,-1/2/dp,Ny)';
[fx,fy]=meshgrid(fx,fy);    % cordinate in spatial frequency domain
AS=fftshift(fft2(fftshift(sig))); % frequency of signal as angular spectrum
% figure,surfl(fx,fy,AS.*conj(AS)),shading interp,colormap(gray); xlabel('fx');ylabel('fy'); title('fft(sig)') %/max(abs(freq(:)))


% ================== NUFFT parameters ================
eps = 10^(-6);  % accuracy of NUFFT
x=reshape(x/(max(abs(x(:))))*pi,Ny*Nx,1);
y=reshape(y/(max(abs(y(:))))*pi,Ny*Nx,1);   % shape 2D data to 1D
%% 0 CV method (T-FFT): based on impluse reponse function 
method0='CV(TFFT)';
sig0=padarray(sig,[Nx/2,Ny/2]);  % zero-padding to double size of source window
[Ny0,Nx0]=size(sig0);
Lx0 = (Nx0)*dp;                  
Ly0 = (Ny0)*dp;
AS0=fftshift(fft2(fftshift(sig0)));
x0 = linspace(-Lx0/2,Lx0/2,Nx0)';              % cordinate in spatial domain
y0 = linspace(Ly0/2,-Ly0/2,Ny0)'; 
[x0,y0]=meshgrid(x0,y0);
r = sqrt(x0.^2+y0.^2+z^2); 
kernel = exp(1i*k*r)./1j/lam./r;  
kernel_FT = fftshift(fft2(fftshift(kernel)));
E0 = ifftshift(ifft2(ifftshift(kernel_FT.*AS0)));
E0  = E0./max(abs(E0(:)));
sft=1;
E0 = E0(Nx0/2-Nx/2+sft:Nx0/2+Nx/2+sft-1,Ny0/2-Ny/2+sft:Ny0/2+Ny/2+sft-1);
figure,imshow(abs(E0),[]);title([method0,' amplitude of diffractive field']),colormap(turbo)
figure,imshow(angle(E0),[]);title([method0,' pahse of diffractive field']),colormap(turbo)
%% I. BL-AS method: [K. Matsushima and T. Shimobaba, Opt. Express17, 19662 (2009)]
% method1='BL-AS';
% tic
% sig1=padarray(sig,[Nx/2,Ny/2]);  % zero-padding to double size of source window
% [Ny1,Nx1]=size(sig1);
% AS1=fftshift(fft2(fftshift(sig1))); % frequency of signal as angular spectrum
% dfx1=1/2/Lx;dfy1=1/2/Ly;
% fx1 = linspace(-1/2/dp,1/2/dp-dfx1,Nx1)';   
% fy1 = linspace(1/2/dp,-1/2/dp+dfy1,Ny1)';
% [fx1,fy1]=meshgrid(fx1,fy1);    % cordinate in spatial frequency domain
% fx_BL = Nx1/2*dp/lam/z;
% fy_BL = Ny1/2*dp/lam/z;    % boundary frequency
% H = exp(1i*k*z*sqrt(1-(fx1*lam).^2-(fy1*lam).^2));
% H = H.*(abs(fx1)<=fx_BL).*(abs(fy1)<=fy_BL);
% FH= AS1.*H;
% E1=ifftshift(ifft2(ifftshift(FH)));
% E1=E1/max(E1(:));
% sft=1;
% E1 = E1(Nx1/2-Nx/2+sft:Nx1/2+Nx/2+sft-1,Ny1/2-Ny/2+sft:Ny1/2+Ny/2+sft-1);
% time1=toc;
% 
% strtime1=strcat(method1,' : ',num2str(time1),'s');
% % figure,imshow(abs(E1),[]);title(method1),colormap(turbo)
% 
% % holo=E1;
% % reconstruct
% % reco1=reco;
% figure,imshow(abs(E1),[]);title(strcat(method1,"reco")),colormap(turbo)
% figure,imshow(angle(E1),[]);title(strcat(method1,"reco")),colormap(turbo)
% % figure,plot(reco1(:,200)),title([method1,'recon']),
%% II. Wide-window (WW-AS) method: [X. Yu, et al., Opt. Lett.37, 4943 (2012)]
% WWAS
%% III. Non-uniform sampling (NuS-AS) method: [Y.-H. Kim,et al., J. Opt.16, 125710 (2014)]
% method3="NuS-AS";
% tic
% fx_BL = Nx*dp/lam/z;
% fy_BL = Ny*dp/lam/z;    % boundary frequency
% Nfx=Nx*2;
% Nfy=Ny*2;
% dfx3=2*fx_BL/(Nfx);
% dfy3=2*fy_BL/(Nfy); % frquency sampling interval
% fxn=linspace(-fx_BL,fx_BL-dfx3,Nfx)';
% fyn=linspace(fy_BL,-fy_BL+dfy3,Nfy)';    % new frequency incoordinates
% [fxn,fyn]=meshgrid(fxn,fyn);
% Diffr_NUFFT2;   % call the function of diffraction progress
% E3=E;
% time3=time;
% 
% strtime3=strcat(method3,':',num2str(time3),'s');
% % figure,imshow(abs(E));title(method3) %
% Nfx3=Nfx;
% Nfy3=Nfy;
%% IV. AS-AS method: [W. Zhang, H. Zhang, and G. Jin, Opt. Lett.45, 4416 (2020)]
method4="AS-AS";
tic;
fx_BL = Nx*dp/lam/z;
fy_BL = Ny*dp/lam/z;    % boundary frequency
Nfx=round(4*fx_BL*Nx*dp);
Nfy=round(4*fy_BL*Ny*dp);  % N_BL  
if Nfx<2 || Nfy<2
    Nfx=2;Nfy=2;  % ensure NUFFT can work
end
dfx4=2*fx_BL/(Nfx);
dfy4=2*fy_BL/(Nfy);   % frquency sampling interval
fxn=linspace(-fx_BL,fx_BL,Nfx)';
fyn=linspace(fy_BL,-fy_BL,Nfy)';    % new frequency coordinates
[fxn,fyn]=meshgrid(fxn,fyn);
Diffr_NUFFT2;   % call the function of diffraction progress
E4=E;   % diffraction field
time4=time;

strtime4=strcat(method4,': ',num2str(time4),'s');
Nfx4=Nfx;
Nfy4=Nfy;
figure,imshow(abs(E4),[]);title(strcat(method4," amplitude of diffractive field")),colormap(turbo)
figure,imshow(angle(E4),[]);title(strcat(method4," phase of diffractive field")),colormap(turbo)
%% V. BE-AS method: [W. Zhang, H. Zhang, and G. Jin, Opt. Lett.45, 1543 (2020)]
method5="BE-AS";
tic
fx_BE=sqrt(Nx/lam/z/2); 
fy_BE=sqrt(Ny/lam/z/2); % boundary frequency
Nfx=Nx*2;
Nfy=Ny*2;   % N_BE
dfx5=2*fx_BE/(Nfx);
dfy5=2*fy_BE/(Nfy);   % frquency sampling interval
fxn=linspace(-fx_BE,fx_BE,Nfx)';
fyn=linspace(fy_BE,-fy_BE,Nfy)';    % new frequency incoordinates
[fxn,fyn]=meshgrid(fxn,fyn);
Diffr_NUFFT2;   % call the function of diffraction progress
E5=E;
time5=time;

strtime5=strcat(method5,': ',num2str(time5),'s');
Nfx5=Nfx;
Nfy5=Nfy;

figure,imshow(abs(E5),[]);title(strcat(method5," amplitude of diffractive field")),colormap(turbo)
figure,imshow(angle(E5),[]);title(strcat(method5," phase of diffractive field")),colormap(turbo)
%% VI. CE-AS method: the proposed method
method6="CE-AS";
%% based on f_BE
AS=AS.*(abs(fx)<=fx_BE).*(abs(fy)<=fy_BE); % one can change fx_BE to fx_BL
ESD=AS.*conj(AS);	% engergy spectrum density
jx=find(fx(1,:)-fx_BE>=0,1,'first');
jy=find(fy(:,1)-fy_BE>=0,1,'last');
Tegy=sum(sum(ESD.*dfx).*dfy);
egy=Tegy;
weight=1;
in=0;
while weight>=eta
    egy=egy-(sum(ESD(jy,Nx+1-jx:jx)*dfx*dfy)*2+sum(ESD(jy:Ny+1-jy,jx)*dfy*dfx)*2-ESD(jy,jx)*dfx*dfy*4);
    jx=jx-1; 
    jy=jy+1; % expand the range to search
    weight=egy/Tegy;
    in=in+1;
end
tic
fx_CE=fx(1,jx);
fy_CE=fy(jy,1);

Nfx=round(fx_CE^2*4*lam*z);
Nfy=round(fy_CE^2*4*lam*z);  % N_SE

dfx6=2*fx_CE/(Nfx);
dfy6=2*fy_CE/(Nfy);   % new frequency sampling interval
fxn=linspace(-fx_CE,fx_CE,Nfx)';
fyn=linspace(fx_CE,-fx_CE,Nfy)';
[fxn,fyn]=meshgrid(fxn,fyn);
Diffr_NUFFT2;
E6=E;
time6=time;

strtime6=strcat(method6,': ',num2str(time6),'s');
Nfx6=Nfx;
Nfy6=Nfy;

figure,imshow(abs(E6),[]);title(strcat(method6," amplitude of diffractive fieldr")),colormap(turbo)
figure,imshow(angle(E6),[]);title(strcat(method6," phase of diffractive field")),colormap(turbo)
%% Display result
E_rsi=E0;
% SNR1=snr(abs(E1),abs(E_rsi)-abs(E1));
% SNR2=snr(abs(E2),abs(E0)-abs(E2));
% SNR3=snr(abs(E3),abs(E_rsi)-abs(E3));
SNR4=snr(abs(E4),abs(E_rsi)-abs(E4));
SNR5=snr(abs(E5),abs(E_rsi)-abs(E5));
SNR6=snr(abs(E6),abs(E_rsi)-abs(E6));

disp(['z=',num2str(zn),'*zc'])
disp('--------------------------------------');
% disp(strcat(method1,': ',num2str(SNR1),' dB'));
% disp(strcat(method2,': ',num2str(SNR2),' dB'));
% disp(strcat(method3,': ',num2str(SNR3),' dB'));
disp(strcat(method4,': ',num2str(SNR4),' dB'));
disp(strcat(method5,': ',num2str(SNR5),' dB'));
disp(strcat(method6,': ',num2str(SNR6),' dB'));

disp('-----------------------------------');
% disp(strtime1);
% disp(strtime2);
% disp(strtime3);
disp(strtime4);
disp(strtime5);
disp(strtime6);
disp('-----------------------------------');
disp(strcat('eta = ',num2str(round(weight,3))));
disp(strcat('accelerate = ',num2str(round(time5/time6,1))));



