% clear all; close all; clc
% clear all; close all
clear all
r=1;

% while 1
    path=['/Volumes/TOSHIBA EXT/Dropbox/life_res/lattice_gases/LBMs/D2Q9/hydro_vent/shear_obst/opts_hyv_' num2str(r) '.txt'];
%     path=['/Users/stuartbartlett/Dropbox/life_res/lattice_gases/LBMs/D2Q9/hydro_vent/shear_obst/opts_hyv_' num2str(r) '.txt'];
%     path=['/Users/stuartbartlett/Dropbox/life_res/lattice_gases/LBMs/D2Q9/thermal/all_mod/opts_sta_' num2str(r) '.txt'];
    path=['/Volumes/TOSHIBA EXT/Dropbox/life_res/lattice_gases/LBMs/D2Q9/thermal/all_mod/carnot/sq_obs/y_sink_1_2/opts_hv_' num2str(r) '.txt'];
%     path=['/Users/stuartbartlett/Dropbox/life_res/lattice_gases/LBMs/D2Q9/thermal/all_mod/carnot/sq_obs/yopts_hv_' num2str(r) '.txt'];
%     path=['/Users/stuartbartlett/Dropbox/life_res/lattice_gases/LBMs/D2Q9/thermal/radiative/opts_' num2str(r) '.txt'];
    
    fid=fopen(path);
    a=fscanf(fid,'%f');
    
    y_s=a(1); x_s=a(2);
    t_end=a(3); t_int=a(4);
    tau_v=a(5); tau_c=a(6);
    b_g0=a(7);
    n_cs=0;
    
    for t=1:t_end/t_int
        for j=1:x_s
            rho(:,j,t)=a(8+((t-1)*(4+n_cs)+0)*y_s*x_s+(j-1)*y_s:7+((t-1)*(4+n_cs)+0)*y_s*x_s+j*y_s);
            eps(:,j,t)=a(8+((t-1)*(4+n_cs)+1)*y_s*x_s+(j-1)*y_s:7+((t-1)*(4+n_cs)+1)*y_s*x_s+j*y_s);
            u_x(:,j,t)=a(8+((t-1)*(4+n_cs)+2)*y_s*x_s+(j-1)*y_s:7+((t-1)*(4+n_cs)+2)*y_s*x_s+j*y_s);
            u_y(:,j,t)=a(8+((t-1)*(4+n_cs)+3)*y_s*x_s+(j-1)*y_s:7+((t-1)*(4+n_cs)+3)*y_s*x_s+j*y_s);
        end
    end
    r=r+1;
% end
clear a

%%

% clear all; close all; clc
% clear all; close all
clear all
r=10;

path=['/Volumes/TOSHIBA EXT/Dropbox/life_res/lattice_gases/LBMs/D2Q9/hydro_vent/shear_obst/opts_hyv_' num2str(r) '.txt'];
% path=['/Users/stuartbartlett/Dropbox/life_res/lattice_gases/LBMs/D2Q9/hydro_vent/shear_obst/opts_hyv_' num2str(r) '.txt'];

fid=fopen(path);
a=fscanf(fid,'%f');

y_s=a(1); x_s=a(2);
t_end=a(3); t_int=a(4);
tau_v=a(5); tau_c=a(6);
b_g0=a(7); 

for j=1:x_s
    rho(:,j)=a(8+(j-1)*y_s:7+j*y_s);
    eps(:,j)=a(8+y_s*x_s+(j-1)*y_s:7+y_s*x_s+j*y_s);
    u_x(:,j)=a(8+2*y_s*x_s+(j-1)*y_s:7+2*y_s*x_s+j*y_s);
    u_y(:,j)=a(8+3*y_s*x_s+(j-1)*y_s:7+3*y_s*x_s+j*y_s);
end
clear a

rho(rho==0)=1;
eps(isnan(eps))=0.5;
close all
figure('units','normalized','outerposition',[0.3 0.03 0.45 0.95])

del_T=eps(end,1)-eps(1,1);
%     Ra=round(b_g0*del_T*H^3/(chi*nu));
%     Nu(t)=1+H*sum(sum(u_y(:,:,t).*eps(:,:,t)/R))/(chi*del_T*o_n);

subplot(4,1,1)
imagesc(rho)
daspect([1 1 1])
text(-0.2*size(rho,2),0.2*size(rho,1),'\rho','fontsize',28)
colorbar

subplot(4,1,2)
imagesc(u_x)
daspect([1 1 1])
text(-0.2*size(rho,2),0.2*size(rho,1),'u','fontsize',28)
colorbar

subplot(4,1,3)
imagesc(u_y)
daspect([1 1 1])
text(-0.2*size(rho,2),0.2*size(rho,1),'w','fontsize',28)
colorbar

subplot(4,1,4)
imagesc(eps)
daspect([1 1 1])
text(-0.2*size(rho,2),0.2*size(rho,1),'\epsilon','fontsize',28)
colorbar

%     text(0.1*size(rho,2),-0.1*size(rho,1),['Ra=' num2str(Ra(t)) ' Nu= ' num2str(Nu(t))],'fontsize',28)
    
    
    
    % t=size(eps,3)-1;
%
% % u_p=u_x(:,:,t);
% % w_p=u_y(:,:,t);
%
% o_set=4;
% 
% for i=o_set+1:y_s-o_set
%     for j=o_set+1:x_s-o_set
%         u_b=reshape(u_x(i-o_set:i+o_set,j-o_set:j+o_set,t),[1 (2*o_set+1)^2]);
%         u_b(u_b==0)=[];        
%         u_p(i-o_set,j-o_set)=mean(u_b);
%         
%         u_b=reshape(u_y(i-o_set:i+o_set,j-o_set:j+o_set,t),[1 (2*o_set+1)^2]);
%         u_b(u_b==0)=[];
%         w_p(i-o_set,j-o_set)=mean(u_b);
%     end
% end
% 
% u_p=[zeros([y_s o_set]) [zeros([o_set x_s-2*o_set]); u_p; zeros([o_set x_s-2*o_set])] zeros([y_s o_set])];
% w_p=[zeros([y_s o_set]) [zeros([o_set x_s-2*o_set]); w_p; zeros([o_set x_s-2*o_set])] zeros([y_s o_set])];
% 
% u_p(isnan(u_p))=0;
% w_p(isnan(w_p))=0;

%%

% o_n=sum(sum(rho(:,:,end-1)~=0));
% rho(rho==0)=rho(1,1,1);
% eps(eps==0)=0.9*(eps(1));
% eps(isnan(eps))=0.9*(eps(1));
% % eps_2=eps;
% % eps_2(eps_2==0)=0.9;
% % % % psi(find(psi==0))=1;
% % 
% % eps_2=(eps_2-min(min(min(eps_2))))/(max(max(max(eps_2)))-min(min(min(eps_2))));
% % u_x_p=(u_x-min(min(min(u_x))))/(max(max(max(u_x)))-min(min(min(u_x))));
% % u_y_p=(u_y-min(min(min(u_y))))/(max(max(max(u_y)))-min(min(min(u_y))));
% 
% % dux_dy=[-u_x(1,:,:)-u_x(2,:,:)/3; 0.5*(u_x(1:end-2,:,:)-u_x(3:end,:,:)); u_x(end-1,:,:)/3+u_x(end,:,:)];
% % duy_dx=0.5*([u_y(:,2:end,:) u_y(:,1,:)]-[u_y(:,end,:) u_y(:,1:end-1,:)]);
% % vort=(dux_dy-duy_dx)/c;

close all
figure('units','normalized','outerposition',[0.3 0.03 0.45 0.95])
vid=0;
u_mag=sqrt(u_x.^2+u_y.^2);

set(gcf,'color','k')
% colormap jet

if vid
    writerObj=VideoWriter(['/Users/affe_prime/Desktop/circ.avi']);
%     writerObj=VideoWriter(['/Users/stuartbartlett/Desktop/sta.avi']);
    writerObj.FrameRate=20;
    writerObj.Quality=100;
    open(writerObj);
end

T_0=0.5*(eps(1,1,1)+eps(end,1,1));
R=1;
c=sqrt(3*R*T_0);
H=y_s*c;
nu=(1/3)*(tau_v-0.5)*c^2;
chi=(2/3)*(tau_c-0.5)*c^2;

% d_fl_t=squeeze(sum(-chi*(eps(1,:,:)-eps(3,:,:)),2)/(2*c^2));
% d_fl_b=squeeze(sum(-chi*(eps(y_s-2,:,:)-eps(y_s,:,:)),2)/(2*c^2));

for t=1:size(eps,3)
%     del_T=eps(end,1,t)/R-eps(1,1,t)/R;
%     Ra(t)=round(b_g0*del_T*H^3/(chi*nu));
%     Nu(t)=1+H*sum(sum(u_y(:,:,t).*eps(:,:,t)/R))/(chi*del_T*o_n);
    
    subplot(4,1,1)
    imagesc(u_mag(:,:,t))
    daspect([1 1 1])
    text(-0.2*size(rho,2),0.4*size(rho,1),'|u+w|','fontsize',28,'color','w')
    
    subplot(4,1,2)
    imagesc(u_x(:,:,t))
    daspect([1 1 1])
    text(-0.2*size(rho,2),0.4*size(rho,1),'u','fontsize',28,'color','w')
    
    subplot(4,1,3)
    imagesc(u_y(:,:,t))
    daspect([1 1 1])
    text(-0.2*size(rho,2),0.4*size(rho,1),'w','fontsize',28,'color','w')
    
    subplot(4,1,4)
    imagesc(eps(:,:,t)/R)
    daspect([1 1 1])
    text(-0.2*size(rho,2),0.4*size(rho,1),'\epsilon','fontsize',28,'color','w')
%     text(0.1*size(rho,2),-0.1*size(rho,1),['Ra=' num2str(Ra(t)) ' Nu= ' num2str(Nu(t))],'fontsize',28)
    
    if vid
        frame=getframe(gcf);
        writeVideo(writerObj,frame);
    end
    
    pause(0.001)
end

if vid
    close(writerObj);
end

%%
close all
x_frac=0.2;
y_frac=0.5;
figure('units','normalized','outerposition',[0.3 0.03 0.45 0.95])
u_mag=sqrt(u_x.^2+u_y.^2);
t=size(rho,3);
subplot(4,2,1)
imagesc(u_mag(:,:,t))
hold
plot([x_frac*x_s x_frac*x_s],[0 y_s],'-r')
plot([0 x_s],[y_frac*y_s y_frac*y_s],'-r')
daspect([1 1 1])
text(0.3*size(rho,2),-0.15*size(rho,1),'sqrt(u^2+w^2)','fontsize',28,'color','k')

subplot(4,2,2)
imagesc(u_mag(:,:,t))
hold
plot([x_frac*x_s x_frac*x_s],[0 y_s],'-r')
plot([0 x_s],[y_frac*y_s y_frac*y_s],'-r')
daspect([1 1 1])
text(0.3*size(rho,2),-0.15*size(rho,1),'sqrt(u^2+w^2)','fontsize',28,'color','k')

subplot(4,2,3)
imagesc(u_x(:,:,t))
hold
plot([x_frac*x_s x_frac*x_s],[0 y_s],'-r')
plot([0 x_s],[y_frac*y_s y_frac*y_s],'-r')
daspect([1 1 1])
text(0.3*size(rho,2),-0.15*size(rho,1),'u','fontsize',28,'color','k')

subplot(4,2,4)
imagesc(abs(u_x(:,:,t)))
hold
plot([x_frac*x_s x_frac*x_s],[0 y_s],'-r')
plot([0 x_s],[y_frac*y_s y_frac*y_s],'-r')
daspect([1 1 1])
text(0.3*size(rho,2),-0.15*size(rho,1),'|u|','fontsize',28,'color','k')

subplot(4,2,5)
imagesc(u_y(:,:,t))
hold
plot([x_frac*x_s x_frac*x_s],[0 y_s],'-r')
plot([0 x_s],[y_frac*y_s y_frac*y_s],'-r')
daspect([1 1 1])
text(0.3*size(rho,2),-0.15*size(rho,1),'w','fontsize',28,'color','k')

subplot(4,2,6)
imagesc(abs(u_y(:,:,t)))
hold
plot([x_frac*x_s x_frac*x_s],[0 y_s],'-r')
plot([0 x_s],[y_frac*y_s y_frac*y_s],'-r')
daspect([1 1 1])
text(0.3*size(rho,2),-0.15*size(rho,1),'|w|','fontsize',28,'color','k')

subplot(4,1,4)
imagesc(eps(:,:,t))
hold
plot([x_frac*x_s x_frac*x_s],[0 y_s],'-r')
plot([(x_frac+0.5)*x_s (x_frac+0.5)*x_s],[0 y_s],'-r')
plot([0 x_s],[y_frac*y_s y_frac*y_s],'-r')
daspect([1 1 1])
text(0.3*size(rho,2),-0.15*size(rho,1),'T','fontsize',28,'color','k')
%     text(0.1*size(rho,2),-0.1*size(rho,1),['Ra=' num2str(Ra(t)) ' Nu= ' num2str(Nu(t))],'fontsize',28)

%%
close all
u_mag=sqrt(u_x.^2+u_y.^2);
figure('units','normalized','outerposition',[0.1 0.2 0.8 0.8])
hold
x_frac=0.2;
plot(0.5*(u_x(140,:,end)+u_x(141,:,end)))
plot(0.5*(u_y(140,:,end)+u_y(141,:,end)))
plot(0.5*(u_mag(140,:,end)+u_mag(141,:,end)))
plot([ceil(x_frac*x_s) ceil(x_frac*x_s)],[-5e-2 5e-2],'-b')
plot([ceil((x_frac+0.25)*x_s) ceil((x_frac+0.25)*x_s)],[-5e-2 5e-2],'-r')
plot([ceil((x_frac+0.5)*x_s) ceil((x_frac+0.5)*x_s)],[-5e-2 5e-2],'-b')
grid on
legend('u','w','|u+w|')

%%

figure
hold
plot(0.5*(u_x(35,:,end)+u_x(36,:,end)))
plot(0.5*(u_y(35,:,end)+u_y(36,:,end)))
plot(0.5*(u_mag(35,:,end)+u_mag(36,:,end)))
plot([ceil(x_frac*x_s) ceil(x_frac*x_s)],[-5e-2 5e-2],'-b')
plot([ceil(0.38*x_s) ceil(0.38*x_s)],[-5e-2 5e-2],'-r')
plot([ceil(0.63*x_s) ceil(0.63*x_s)],[-5e-2 5e-2],'-b')
grid on
legend('u','w','|u+w|')

%%

close all

figure('units','normalized','outerposition',[0.3 0.05 0.38 0.9])
set(gcf,'color','w')
colormap jet

subplot(4,1,1)
imagesc(P(:,:,t)*3/c^2)
% imagesc(P(4:end-3,4:end-3,t))
daspect([1 1 1])
text(-250,200,'P[bar]','fontsize',28)
set(gca,'fontsize',18)
colorbar


% stream particles animation 
% 
% diagram showing velocity magnitude to show dead zone

% streamslice(repmat(1:x_s,[240 1]),repmat((1:240)',[1 x_s]),u_x(1:240,:,end-1),u_y(1:240,:,end-1),1)

% streamslice(repmat(1:x_s-6,[y_s-6 1]),repmat((1:y_s-6)',[1 x_s-6]),u_p,w_p,1)
streamslice(repmat(1:x_s,[y_s 1]),repmat((1:y_s)',[1 x_s]),u_p,w_p,1.3)


% 
% 
% sx=ones([20 1]);
% sy=linspace(0.1*y_s,0.6*y_s,20);
% 
% streamline(stream2(repmat(1:x_s,[y_s 1]),repmat((1:y_s)',[1 x_s]),u_x(:,:,end-1),u_y(:,:,end-1),sx,sy));
% 
% sx=repmat(600,[1 10]);
% 
% sx=linspace(500,600,10)
% sy=linspace(165,240,10);
% 
% 
% for i=1:10
%     sx=[sx linspace(500+55*i,600+55*i,10)];
%     sy=[sy linspace(165,240,10)];
% end
% 
% % sy=linspace(165,240,10);
% 
% streamline(stream2(repmat(1:x_s,[y_s 1]),repmat((1:y_s)',[1 x_s]),u_x(:,:,end-1),u_y(:,:,end-1),sx,sy));
% 
% daspect([1 1 1])
% set(gca,'ydir','reverse')


subplot(4,1,2)
imagesc(u_x(:,:,t)/c)
daspect([1 1 1])
text(-140,200,'u','fontsize',28)
set(gca,'fontsize',18)
colorbar

subplot(4,1,3)
imagesc(u_y(:,:,t)/c)
daspect([1 1 1])
text(-140,200,'w','fontsize',28)
set(gca,'fontsize',18)
colorbar

subplot(4,1,4)
imagesc(eps(:,:,t)-273.15)
daspect([1 1 1])
text(-210,200,'T[C]','fontsize',28)
set(gca,'fontsize',18)
colorbar



%%

sx=ones([20 1]);
sy=linspace(0.1*y_s,0.6*y_s,20);


figure
imagesc(P(:,:,end-1))
streamline(stream2(repmat(1:x_s,[y_s 1]),repmat((1:y_s)',[1 x_s]),u_x(:,:,end-1),u_y(:,:,end-1),sx,sy));
daspect([1 1 1])
set(gca,'ydir','reverse')



