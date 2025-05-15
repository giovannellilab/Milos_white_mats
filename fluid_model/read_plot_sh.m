clear all; close all; clc
% clear all; close all
clear all

path='/Users/stuartbartlett/Dropbox/life_res/lattice_gases/LBMs/D2Q9/hydro_vent/shear_obst/opts_hy_v.txt';

fid=fopen(path);
a=fscanf(fid,'%f');

y_s=a(1); x_s=a(2);
t_end=a(3); t_int=a(4);
n_cs=a(5);
k=6;

for t=1:t_end/t_int
    for j=1:x_s
        for i=1:x_s
            eps(i,j,t)=a(k+(3+n_cs)*((t-1)*y_s*x_s+((j-1)*y_s+i)));
            u_x(i,j,t)=a(k+(3+n_cs)*((t-1)*y_s*x_s+((j-1)*y_s+i))+1);
            u_y(i,j,t)=a(k+(3+n_cs)*((t-1)*y_s*x_s+((j-1)*y_s+i))+2);
            for k=1:n_cs
                psi(i,j,k,t)=a(k+(3+n_cs)*((t-1)*y_s*x_s+((j-1)*y_s+i))+2+k);
            end
        end
    end
end
clear a

%%
close all
figure
colormap jet
for t=1:t_end/t_int
    subplot(2,1,1)
    imagesc(u_x(:,:,t)/U)
    daspect([1 1 1])
    
    subplot(2,1,2)
    imagesc(u_y(:,:,t)/U)
    daspect([1 1 1])
    
    text(470,-75,['t=' num2str(t)])
    pause(0.02)
end


