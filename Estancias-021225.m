% Definir o declarar varibles
Ar1='A3.nc';
Run1=' ';
Var1='U';
Var2='T';
Var3='W';
time=43;
dx=10;
puntosx=2001;
puntosz=99;
ct=28
% Definir ejex
dis=(puntosx-1)*dx/2;
ejex1=((-dis):dx:(dis))/1000;
clear dis 
ejex2=ejex1;
ejex1(length(ejex1))=[];

gravedad=1.352;
timbebar=1:1:time;
PHB1=squeeze(ncread(Ar1,'PHB'));
PH1=squeeze(ncread(Ar1,'PH'));
gz=PHB1+PH1;
clear PHB PH

gz2d=mean(gz,3,"omitmissing");
ejez1=mean(gz2d)./ (gravedad);
for i= 1:1:puntosz
    ejez(i)=(ejez1(i)+ejez1(i+1))/2;
end
clear gz2d gz

un=squeeze(ncread("A3.nc",'U'));
t=squeeze(ncread("A3.nc",'T'))+100;

for i = 1:1:puntosx-1
    u(i,:,:)=(un(i,:,:)+un(i+1,:,:))/2;
end

figure;
contourf(ejex1,ejez,t(:,:,ct)',500,'EdgeColor','none');
colormap(turbo), colorbar;
c=colorbar;
mincolor=min(mean(t(:,:,ct)));
maxcolor=max(mean(t(:,:,ct)));
clim([mincolor,maxcolor])
ylabel('Height (m)')
xlabel('Distance (Km)')
hold on 

contour(ejex1,ejez,u(:,:,ct)', 'EdgeColor','black')















