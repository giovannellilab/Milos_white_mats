clear all; close all; clc

b_m=1.22;
b_std=0.37;
b_sk=2.4;
b_kt=15;
% Measured sand grain distribution
gr_sd=[0.36 0.009; 0.55 0.019; 0.74 0.031; 0.93 0.288; 1.13 0.244; 1.32 0.153; 1.51 0.072; 1.71 0.028; 1.90 0.022; 2.09 0.013; 2.29 0.009; 2.48 0.003; 2.67 0.003; 2.87 0; 3.06 0.006];
% Normalise
gr_sd(:,2)=gr_sd(:,2)/(sum(gr_sd(:,2))*mean(diff(gr_sd(:,1))));
% % Optimise distribution parameters
% err=1e7;
% for gr_m=1.21:0.001:1.23
%     for gr_std=0.36:0.001:0.38
%         for skew=2.3:0.01:2.5
%             for kurt=14:0.1:16
%                 gr_samp=pearsrnd(gr_m,gr_std,skew,kurt,[2e5 1]);
%                 [counts ctrs]=hist(gr_samp,gr_sd(:,1));
%                 counts=counts/(sum(counts)*mean(diff(ctrs)));
%                 err_n=sum(abs(counts-gr_sd(:,2)'));
%                 if err_n<err
%                     err=err_n;
%                     b_m=gr_m;
%                     b_std=gr_std;
%                     b_sk=skew;
%                     b_kt=kurt;
%                 end
%             end
%         end
%     end
% end

% Plot optimal distribution
gr_samp=pearsrnd(b_m,b_std,b_sk,b_kt,[1e6 1]);
[counts ctrs]=hist(gr_samp,linspace(0,4,100));
counts=counts/(sum(counts)*mean(diff(ctrs)));

% close all
% figure
% hold
% plot(ctrs,counts,'b-','linewidth',1.5)
% plot(gr_sd(:,1),gr_sd(:,2),'go-','markersize',8,'linewidth',0.5)
% xlabel('Grain size x [mm]')
% ylabel('\phi(x)')
% xlim([0 3.3])
% set(gca,'fontsize',14)
% set(gca,'xtick',0:0.5:3)

sc=20;              % sc=grid cells per mm
x_s=10080;          % 1e4 points=0.5m, grid cells are 0.05mm across
psty=0.55;          % Porosity

% sc=2;               % sc=grid cells per mm
% x_s=1008;           % 1e3 points=0.5m, grid cells are 0.5mm across
% psty=0.7;          % Porosity

y_s=25*10*sc;       % 25cm total depth
tri_w=8*10*sc;      % 8cm wide
tri_h=5*10*sc;      % 5cm tall
g_olap=0.2;         % Grain overlap
ch_dep=15*10*sc;    % Channel depth
t_st=10*10*sc;      % x position at which ripple is placed
% Grid pad sizes
x_p=round(0.1*x_s);
y_p=round(0.1*y_s);
% Initialise grid for indexing porous region
p_grid=false([y_s+2*y_p x_s+2*x_p]);
% Index lower porous region
% p_grid(y_p+ch_dep+1:y_s+y_p-round(0.02*y_s),x_p+1:x_s+x_p)=1;
p_grid(y_p+ch_dep+1:y_s+y_p,x_p+1:x_s+x_p)=1;
% Define triangle region
m_sec=false([tri_h tri_w]);
for i=tri_h:-1:1
    m_sec(i,1+round(0.5*tri_w*(1-i/tri_h)):round(0.5*tri_w*(1+i/tri_h)))=1;
end
% Add triangle region
p_grid(y_p+ch_dep-tri_h+1:y_p+ch_dep,x_p+t_st+1:x_p+t_st+tri_w)=m_sec;
% Display porous media index grid
close all
figure
imagesc(p_grid(y_p+1:y_p+y_s,x_p+1:x_p+x_s))
daspect([1 1 1])
% Collect indices for points where grains can be placed
[ms ns]=find(p_grid);
inds=find(p_grid);
p_grid=false([y_s+2*y_p x_s+2*x_p]);

%%
close all
% Construct porous medium
while 1
    % Generate polygon side number
    gr_s=round(randn+6);
    % Eliminate if less than 4 or greater than 8
    gr_s=4*(gr_s<4)+gr_s*(gr_s<=8 && gr_s>=4)+8*(gr_s>8);
    % Generate grain diameter
    gr_d=pearsrnd(b_m,b_std,b_sk,b_kt,1);
    if gr_d>=4
        continue
    end
    % Generate lengths of vertex points. Add normally distributed, 5% noise to lengths
    lens=0.5*gr_d*(1+0.05*randn([gr_s 1]));
    lens(lens<0)=b_m;
    % Generate angles of vertices
    angs=sort(linspace(0,2*pi*(1-1/gr_s),gr_s)'+0.2*randn([gr_s 1]));
    % Calculate coordinates of vertices
    coords=sc*[lens.*cos(angs) lens.*sin(angs)];
    % Close the loop of points
    coords=[coords; coords(1,:)];
    % Rotate grain by random angle
    ang=2*pi*rand;
    coords=coords*[cos(ang) -sin(ang); sin(ang) cos(ang)];
    % Find bounding box for grain and create BW mask
    mx_c=max(coords);
    mn_c=min(coords);
    bx_m=ceil(mx_c(2)-mn_c(2));
    bx_n=ceil(mx_c(1)-mn_c(1));
    ext=poly2mask(coords(:,1)-mn_c(1)+1,coords(:,2)-mn_c(2)+1,bx_m,bx_n);
    % Generate random position to place grain
    p_i=round((length(ms)-1)*rand)+1;
    c_p=[ns(p_i) ms(p_i)];
    % Check if brain will be placed near the bottom of the grid
    if c_p(2)+bx_m-1>round(0.99*y_s)+y_p
        continue
    end
    bx_c=p_grid(c_p(2):c_p(2)+bx_m-1,c_p(1):c_p(1)+bx_n-1)+ext;
    % Check overlap and place if overlap is below threshold
    if (sum(sum(bx_c==2))/sum(sum(ext)))<g_olap
        p_grid(c_p(2):c_p(2)+bx_m-1,c_p(1):c_p(1)+bx_n-1)=p_grid(c_p(2):c_p(2)+bx_m-1,c_p(1):c_p(1)+bx_n-1) | ext;
        if c_p(1)+bx_n-1>x_p+x_s
            p_grid(c_p(2):c_p(2)+bx_m-1,c_p(1)-x_s:c_p(1)-x_s+bx_n-1)=p_grid(c_p(2):c_p(2)+bx_m-1,c_p(1)-x_s:c_p(1)-x_s+bx_n-1) | ext;
        elseif c_p(1)<x_p+0.02*x_s
            p_grid(c_p(2):c_p(2)+bx_m-1,c_p(1)+x_s:c_p(1)+x_s+bx_n-1)=p_grid(c_p(2):c_p(2)+bx_m-1,c_p(1)+x_s:c_p(1)+x_s+bx_n-1) | ext;
        end
    end
    % Check if porosity is still below threshold
    if sum(p_grid(inds))/length(inds)>=(1-psty)
        break
    end
end

%%
% Remove padding
p_grid=p_grid(y_p+1:y_p+y_s,x_p+1:x_p+x_s);
% Fill holes
p_grid=imfill(p_grid,'holes');
p_grid=-1*(p_grid-1);
% Check porosity
fin_psty=sum(sum(p_grid(ch_dep+1:end,:)))/numel(p_grid(ch_dep+1:end,:))
% Display final grid
figure('units','normalized','outerposition',[0.1 0.2 0.85 0.7])
imagesc(p_grid)
daspect([1 1 1])
% Save grids
save('/Volumes/TOSHIBA EXT/Dropbox/life_res/lattice_gases/LBMs/D2Q9/hydro_vent/shear_obst/p_grid.mat')
% save('/Users/stuartbartlett/Dropbox/life_res/lattice_gases/LBMs/D2Q9/hydro_vent/shear_obst/p_grid.mat')

%%

% p_grid=ones([200 400]);
path='/Volumes/TOSHIBA EXT/Dropbox/life_res/lattice_gases/LBMs/D2Q9/hydro_vent/shear_obst/obst_lrg.txt';
% path='/Users/stuartbartlett/Dropbox/life_res/lattice_gases/LBMs/D2Q9/hydro_vent/shear_obst/obst_lrg.txt';
fid=fopen(path,'w');
fprintf(fid,'%d %d ',size(p_grid,1),size(p_grid,2));
for j=1:size(p_grid,2)
    for i=1:size(p_grid,1)
        fprintf(fid,'%d ',p_grid(i,j));
    end
end
fclose(fid)

%%

close all; clc
figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8])
% figure('units','pixels','outerposition',[200 150 2000 2000*size(p_grid,2)/size(p_grid,1)])
image(64*(1-p_grid))
axis off
hold
daspect([1 1 1])

x_b=[5000 5600];
y_b=[3500 3800];
sc=5;
set(gca,'fontsize',20)
% Small box
plot(x_b([1 2 2 1 1]),y_b([1 1 2 2 1]),'w-','linewidth',3)

x_b2=[0.62*size(p_grid,2) 0.62*size(p_grid,2)+sc*(x_b(2)-x_b(1))];
y_b2=[0.5*size(p_grid,1) 0.5*size(p_grid,1)-sc*(y_b(2)-y_b(1))];
% Large box
plot(x_b2([1 1 2 2 1]),y_b2([2 1 1 2 2]),'w-','linewidth',4)
% Left line
plot([x_b(1) x_b2(1)],[y_b(1) y_b2(2)],'w-','linewidth',3)
% Right line
plot([x_b(2) x_b2(2)],[y_b(2) y_b2(1)],'w-','linewidth',3)
% Centre line
plot([x_b(1) x_b2(1)],[y_b(2) y_b2(1)],'w-','linewidth',3)
% Centre line
% plot([x_b(2) x_b2(2)],[y_b(1) y_b2(2)],'w-','linewidth',3)

fa=get(gca,'Position')
% % axes('Position',[fa(1)+0.62*fa(3)-0.08 fa(2)+0.5*fa(4)-0.005 1.6*fa(3)*sc*(x_b(2)-x_b(1))/size(p_grid,2) 0.99*fa(4)*sc*(y_b(2)-y_b(1))/size(p_grid,1)])
axes('position',[fa(1)+0.62*fa(3)+0.001 fa(2)+0.5*fa(4)-0.011 0.99*fa(3)*sc*(x_b(2)-x_b(1))/size(p_grid,2) 0.99*fa(4)*sc*(y_b(2)-y_b(1))/size(p_grid,1)])
image((64*(1-p_grid(y_b(1):y_b(2),x_b(1):x_b(2)))))
axis off
daspect([1 1 1])
print('-depsc2','/Volumes/TOSHIBA EXT/Dropbox/life_res/lattice_gases/LBMs/D2Q9/hydro_vent/shear_obst/sb.eps')
% % print('-depsc2','/Users/stuartbartlett/Dropbox/life_res/lattice_gases/LBMs/D2Q9/hydro_vent/shear_obst/sb.eps')






%%


for i=1:1e6
    % Generate polygon side number
    gr_s=round(randn+6);
    % Eliminate if less than 4 or greater than 8
    gr_s=4*(gr_s<4)+gr_s*(gr_s<=8 && gr_s>=4)+8*(gr_s>8);
    % Generate grain diameter
    gr_d=pearsrnd(b_m,b_std,b_sk,b_kt,1);
    % Generate lengths of vertex points. Add normally distributed, 5% noise to lengths
    lens=0.5*gr_d*(1+0.05*randn([gr_s 1]));
    lens(lens<0)=b_m;
    % Generate angles of vertices
    angs=sort(linspace(0,2*pi*(1-1/gr_s),gr_s)'+0.2*randn([gr_s 1]));
    
    % Calculate coordinates of vertices
    coords=sc*[lens.*cos(angs) lens.*sin(angs)];
    % Close the loop of points
    coords=[coords; coords(1,:)];
    % Rotate grain by random angle
    ang=2*pi*rand;
    coords=coords*[cos(ang) -sin(ang); sin(ang) cos(ang)];
%     plot(0.1*coords(:,1)+2*mod(i,5)+1,0.1*coords(:,2)+2*ceil(i/5)-1,'-','linewidth',2)
      
    mx_c=max(coords);
    mn_c=min(coords);
    bx_m=ceil(mx_c(2)-mn_c(2));
    bx_n=ceil(mx_c(1)-mn_c(1));
    ext=poly2mask(coords(:,1)-mn_c(1)+1,coords(:,2)-mn_c(2)+1,bx_m,bx_n);
    
    %     close all
    %     figure
    %     imagesc(ext)
    %     hold
    %     plot(coords(:,1)-mn_c(1)+1,coords(:,2)-mn_c(2)+1,'ro-')
    %     daspect([1 1 1])
    
    c_p=[round(x_s*(0.1+rand)) round(y_s*(0.6+0.5*rand))];
    bx_c=p_grid(c_p(2):c_p(2)+bx_m-1,c_p(1):c_p(1)+bx_n-1)+ext;
    if (sum(sum(bx_c==2))/sum(sum(ext)))<0.1
        p_grid(c_p(2):c_p(2)+bx_m-1,c_p(1):c_p(1)+bx_n-1)=p_grid(c_p(2):c_p(2)+bx_m-1,c_p(1):c_p(1)+bx_n-1) | ext;
    end
    s_grid=p_grid(0.6*y_s:1.1*y_s,0.1*x_s:1.1*x_s);
    if sum(sum(s_grid))/numel(s_grid)>=(1-psty)
        break
    end
end

s_grid=p_grid(0.1*y_s:1.1*y_s,0.1*x_s:1.1*x_s);
s_grid=imfill(s_grid,'holes');
1-sum(sum(s_grid(0.5*y_s:end,:)))/(0.5*numel(s_grid))

close all
figure('units','normalized','outerposition',[0.1 0.05 0.6 0.9])
imagesc(s_grid)
daspect([1 1 1])

save('b_grid.mat')

%%

% Set triangular ripple size
tri_w=round(0.3*x_s);
tri_h=round(0.1*x_s);

m_sec=false([round(1.2*tri_h) round(1.2*tri_w)]);
for i=tri_h:-1:1
    m_sec(i+round(0.1*tri_h),1+round(0.5*tri_w*(1-i/tri_h)+0.1*tri_w):round(0.5*tri_w*(1+i/tri_h)+0.1*tri_w))=true;
end

close all
figure
imagesc(m_sec)
daspect([1 1 1])

[ms ns]=find(m_sec);
inds=find(m_sec);
m_sec=false([round(1.2*tri_h) round(1.2*tri_w)]);

%%

for i=1:1e5
    % Generate polygon side number
    gr_s=round(randn+6);
    % Eliminate if less than 4 or greater than 8
    gr_s=4*(gr_s<4)+gr_s*(gr_s<=8 && gr_s>=4)+8*(gr_s>8);
    % Generate grain diameter
    gr_d=pearsrnd(b_m,b_std,b_sk,b_kt,1);
    % Generate lengths of vertex points. Add normally distributed, 5% noise to lengths
    lens=0.5*gr_d*(1+0.05*randn([gr_s 1]));
    lens(lens<0)=b_m;
    % Generate angles of vertices
    angs=sort(linspace(0,2*pi*(1-1/gr_s),gr_s)'+0.2*randn([gr_s 1]));
    
    % Calculate coordinates of vertices
    coords=sc*[lens.*cos(angs) lens.*sin(angs)];
    % Close the loop of points
    coords=[coords; coords(1,:)];
    % Rotate grain by random angle
    ang=2*pi*rand;
    coords=coords*[cos(ang) -sin(ang); sin(ang) cos(ang)];
      
    mx_c=max(coords);
    mn_c=min(coords);
    bx_m=ceil(mx_c(2)-mn_c(2));
    bx_n=ceil(mx_c(1)-mn_c(1));
    ext=poly2mask(coords(:,1)-mn_c(1)+1,coords(:,2)-mn_c(2)+1,bx_m,bx_n);
    
    p_i=round((length(ms)-1)*rand)+1;
    c_p=[ns(p_i) ms(p_i)];
    
    bx_c=m_sec(c_p(2):c_p(2)+bx_m-1,c_p(1):c_p(1)+bx_n-1)+ext;
    
    if (sum(sum(bx_c==2))/sum(sum(ext)))<0.1
        m_sec(c_p(2):c_p(2)+bx_m-1,c_p(1):c_p(1)+bx_n-1)=m_sec(c_p(2):c_p(2)+bx_m-1,c_p(1):c_p(1)+bx_n-1) | ext;
    end
    
    if sum(m_sec(inds))/length(inds)>=(1-psty)
        break
    end
end

m_sec=imfill(m_sec,'holes');
1-sum(m_sec(inds))/length(inds)

close all
figure
imagesc(m_sec(0.1*tri_h:1.1*tri_h,0.1*tri_w:1.1*tri_w))
daspect([1 1 1])

%%



p_grid(0.6*y_s-tri_h+1:0.6*y_s,0.6*x_s-0.5*tri_w+1:0.6*x_s+0.5*tri_w)=



%%

pgrid=[zeros([round(ch_frac*y_s)-tri_h x_s]); m_sec; ones([round((1-ch_frac)*y_s) x_s])];

figure('units','normalized','outerposition',[0.2 0.2 0.7 0.7])
imagesc(pgrid)
daspect([1 1 1])

gr_m=1.2;
gr_std=0.26;
skew=0.7;
kurt=3;

% gr_sd=[0.36   0.9; 0.55   1.9; 0.74   3.1; 0.93   28.8; 1.13   24.4; 1.32   15.3; 1.51   7.2; 1.71   2.8; 1.90   2.2; 2.09   1.3; 2.29   0.9; 2.48   0.3; 2.67   0.3; 2.87   0.0; 3.06   0.6];
% gr_samp=gr_std*randn([1e5 1])+gr_m;
% gr_samp=pearsrnd(gr_m,gr_std,skew,kurt,[1e5 1]);
%
% [counts ctrs]=hist(gr_samp,linspace(min(gr_sd(:,1)),max(gr_sd(:,1)),1e2));
%
% plot(ctrs,counts/max(counts),'bo-')
% hold
% plot(gr_sd(:,1),gr_sd(:,2)/max(gr_sd(:,2)),'go-')

% close all
% figure
% hold
% daspect([1 1 1])
% grid on

posty=0.4;
x_p=round(0.1*x_s);
y_p=round(0.1*y_s);
grd=zeros([2*y_p+y_s 2*x_p+x_s]);

%%

clc
s_grid=zeros([y_s x_s]);
for i=1:1e5
    gr_s=round(0.9*randn+6);
    gr_s=4*(gr_s<4)+gr_s*(gr_s<=8 && gr_s>=4)+8*(gr_s>8);
    a_gr_s(i)=gr_s;
    %     gr_d=gr_std*randn+gr_m;
    gr_d=pearsrnd(gr_m,gr_std,skew,kurt,1)
    angs=sort(linspace(0,2*pi*(1-1/gr_s),gr_s)'+0.09*randn([gr_s 1]));
    lens=0.02*randn([gr_s 1])+0.5*gr_d;
    lens(lens<0)=0.01;
    coords=[lens.*cos(angs) lens.*sin(angs)];
    coords=[coords; coords(1,:)];
    %     plot(coords(:,1)+mod(i,10)+1,coords(:,2)+ceil(i/10),'-')
    ang=2*pi*rand;
    
    coords=coords*[cos(ang) -sin(ang); sin(ang) cos(ang)];
    coords=coords-repmat(mean(coords),[size(coords,1) 1])+repmat([(x_s-1)*rand (y_s-1)*rand],[size(coords,1) 1]);
    
    c_m=ceil(mean(coords,1));
    
    %     if grid(c_m(2),c_m(1))
    ext=poly2mask(coords(:,1)+0.1*x_s,coords(:,2)+0.1*y_s,size(grd,1),size(grd,2));
    g2=grd+ext;
    g2(g2>1)=1;
    if sum(sum(g2-grd))>=round(0.9*sum(sum(ext)))
        grd=g2;
        grd(grd>1)=1;
    end
    %     end
    
    s_grd=grd((y_p+ch_frac*y_s)+1:(y_p+y_s),x_p+1:(x_p+x_s));
    
    if sum(sum(s_grd))/numel(s_grd)>=(1-posty)
        break
    end
end

1-sum(sum(s_grd))/numel(s_grd)

close all
figure('units','normalized','outerposition',[0.1 0.05 0.6 0.9])
imagesc(grd(y_p+1:y_p+y_s,x_p+1:x_p+x_s))
daspect([1 1 1])

%%

pgrid=grd(fac*y_p+1:fac*(y_p+y_s),fac*x_p+1:fac*(x_p+x_s));
pgrid=pgrid(1:2496,1:9984);
pgrid(end,:)=0;
pgrid(:,[1 end])=0;

path='/Volumes/TOSHIBA EXT/Dropbox/life_res/lattice_gases/LBMs/D2Q9/hydro_vent/shear_obst/obst.txt';
fid=fopen(path,'w');
for j=1:size(pgrid,2)
    for i=1:size(pgrid,1)
        fprintf(fid,'%d ',pgrid(i,j));
    end
end
fclose(fid)

%%

clear all; close all; clc

x_s=400;
y_s=x_s/2;

ch_frac=0.7;
tri_w=round(0.075*x_s);
tri_h=round(0.075*x_s);
tri_st=round(0.2*x_s);
posty=0.55;

m_sec=rand([tri_h x_s])>posty;
m_sec(:,1:tri_st)=0;
m_sec(:,tri_st+2*tri_w:end)=0;
for i=1:tri_h
    m_sec(i,tri_st:tri_st+tri_w*(tri_h-i)/tri_h)=0;
    m_sec(i,tri_st+tri_w+tri_w*i/tri_h:tri_st+2*tri_w)=0;
end

pgrid=[zeros([round(ch_frac*y_s)-tri_h x_s]); m_sec; rand([round((1-ch_frac)*y_s) x_s])>posty];

h=fspecial('gaussian',[3 3],1);
pgrid=imfill(im2bw(imfilter(pgrid,h),0.5),'holes');

pgrid([1:2 end-1:end],:)=0;
pgrid(:,[1:2 end-1:end])=0;
%
% % for j=1:20:x_s
% % grid(:,j:j+15)=0;
% % end
% grid(:,x_s/2:end)=0;

pgrid=-1*(pgrid-1);

% grid(grid==0)=1;

% grid=ones([y_s x_s]);
% grid(round(2*y_s/6):round(4*y_s/6),round(3*x_s/8):round(5*x_s/8))=0;

1-sum(sum(pgrid(round(ch_frac*y_s):end,:)))/(ch_frac*numel(pgrid))

imagesc(pgrid)
daspect([1 1 1])


path='/Volumes/TOSHIBA EXT/Dropbox/life_res/lattice_gases/LBMs/D2Q9/hydro_vent/shear_obst/obst.txt';
fid=fopen(path,'w');
for j=1:size(pgrid,2)
    for i=1:size(pgrid,1)
        fprintf(fid,'%d ',pgrid(i,j));
    end
end
fclose(fid)









%%

x_s=480;
y_s=240;

rng(28)
ch_frac=round(0.2*y_s);
posty=0.7;
tri_w=round(0.2*x_s);
tri_h=round(0.65*x_s);

m_sec=rand([tri_h x_s])>posty;
for i=1:tri_h
    m_sec(i,1:0.5*x_s-round(0.5*tri_w*i/tri_h))=0;
    m_sec(i,0.5*x_s+round(0.5*tri_w*i/tri_h):end)=0;
end

pgrid=[];
pgrid=[zeros([y_s-ch_frac-tri_h x_s]); m_sec; rand([ch_frac x_s])>posty];
pgrid=-1*(pgrid-1);
% grid(1,:)=1;
% grid(y_s,:)=1;

% grid=ones([y_s x_s]);
% grid(10:150,10:150)=0;
% grid(:,80:100)=1;
% grid(50:75,:)=1;

% grid(:,[175:215 320:400])=1;
pgrid(end-5:end,:)=1;
pgrid(:,1:2)=1;
pgrid(:,end-1:end)=1;
% grid(:,round(0.492*x_s):round(0.508*x_s))=1;

% for i=1:10
%     y_p=round(0.8*ch_frac*rand);
%     l=round(0.7*x_s*rand);
%     x_p=round(0.9*(x_s-l)*rand);
%     th=round(6*rand);
%
%     grid(y_s-ch_frac+y_p:y_s-ch_frac+y_p+th,x_p:x_p+l)=1;
% end


% grid(end,:)=1;

% grid=ones([y_s x_s]);
close all
imagesc(pgrid)
daspect([1 1 1])

%%


