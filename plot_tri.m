function [ triangle] = plot_tri(V,Lx_p,Ly_p,Nx_p, Ny_p)
x=linspace(-Lx_p/2,Lx_p/2,Nx_p);
y=linspace(Ly_p/2,-Ly_p/2,Ny_p);
[x,y]=meshgrid(x,y);
Ap=V(:,1);
Bp=V(:,2);
Cp=V(:,3);
OAOB=(Ap(1)-x).*(Bp(2)-y)-(Bp(1)-x).*(Ap(2)-y);
OBOC=(Bp(1)-x).*(Cp(2)-y)-(Cp(1)-x).*(Bp(2)-y);
OCOA=(Cp(1)-x).*(Ap(2)-y)-(Ap(1)-x).*(Cp(2)-y);
lgc_AB=OAOB>=0;
lgc_BC=OBOC>=0;
lgc_CA=OCOA>=0;
lgc=lgc_AB+lgc_BC+lgc_CA;
triangle=(lgc==0)+(lgc==3);
end

