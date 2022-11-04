% This function determines the plume radius, peak backscatter intensity,
% and the coordinates of plume center on each horizontal cross-section of 
% gridded COVIS Imaging data using the method of region growing. 

function center = get_plume_centerline_properties_fun_regiongrow_nlinfit(covis,xs,ys,Ro,Ri,zmin,zmax,diagnosis)
% Input:
% covis: Matlab structure that includes gridded Imaging data
% xs: x-coordinate of sourve vent (e.g., Inferno) (m)
% ys: y-coordinate of source vent (m)
% Ri: outer search radius (m)
% Ro: inner search radius (m)
% diagnosis: toggle for turning on and off the diagnostic plots
% Output:
% center: Matlab structure that includes estimated plume properties 



% Example
% xs = -13.27;
% ys = -0.77;
% Ro = 8;
% Ri = 4;
% zmin = 5;
% zmax = 14;
% diagnosis = 1;

global fac
fac = 1;

% extract gridded backscatter intensity and coordinates
xg = covis.grid{3}.x;
yg = covis.grid{3}.y;
zg = covis.grid{3}.z;
xgg = squeeze(xg(:,:,1));
ygg = squeeze(yg(:,:,1));
zgg = squeeze(zg(1,1,:));
dx = median(diff(xgg(1,:)));
dy = median(diff(ygg(:,1)));
vg_a = covis.grid{3}.Ia;
vg_d = covis.grid{3}.Id;
dvg = vg_d./vg_a;
vg = vg_d;
vg(10*log10(dvg)<3) = 10^-9;
vg = smooth3(vg,'gaussian',7,0.3);



% domain of interest
xrange = [xs-Ro,xs+Ro];
yrange = [ys-Ro,ys+Ro];

% main loop
z_vent = 4;
[~,i0] = min(abs(zgg-zmin));
[~,iN] = min(abs(zgg-zmax));
vg_thresh = exp(-8);
x0 = xs;
y0 = ys;

count = 0;
x0_out = zeros(1,iN-i0+1);
y0_out = zeros(1,iN-i0+1);
z0_out = zeros(1,iN-i0+1);
area_out = zeros(1,iN-i0+1);
b_out = zeros(1,iN-i0+1);
I_out = zeros(1,iN-i0+1);
MSE_out = zeros(1,iN-i0+1);
X0n0 = x0;
Y0n0 = y0;
for iz=i0:iN
    count = count+1;
    zlevel = zgg(iz);
    ix=find(xgg(1,:)>=xrange(1)&xgg(1,:)<=xrange(2));
    iy=find(ygg(:,1)>=yrange(1)&ygg(:,1)<=yrange(2));
    % calculate the horizontal index
    X1 = xgg(iy,ix);
    Y1 = ygg(iy,ix);
    % backscatter intensity within the domain of interest
    vg1=vg(iy,ix,iz); % volume scattering coefficient (m^-1)  
    d = sqrt((X1(:)-X0n0).^2+(Y1(:)-Y0n0).^2);
    vg1_search = vg1;
    vg1_search(d>Ri) = 0;
    vg1_base = max(vg1_search(:));
    vg1_search = vg1_search/vg1_base;
    ww = zeros(size(vg1_search));
    ww(vg1_search>vg_thresh)=1;
    cc = bwconncomp(ww,4);
    object_prop = regionprops(cc);
    if isempty(object_prop) || all(ww(:)==0)
        x0_out(count) = X0n;
        y0_out(count) = Y0n;
        z0_out(count) = zlevel;
        area_out(count) = nan;
        I_out(count) = nan;
        b_out(count) = nan;
        continue
    end
    object_area = [object_prop.Area];
    [area_max,i_area]=max(object_area);
    indc = object_prop(i_area).Centroid;
    ind = cc.PixelIdxList;
    indp = ind{i_area};
    mask = zeros(size(vg1));
    mask(indp) = 1;
    vg1_mask = vg1_search.*mask;
    
    % determine centerline properties by finding the best fit of a 2D Gaussian profile
    options = statset('RobustWgtFun','fair');
    ind_fit = find(vg1_mask~=0);
    An0 = quantile(vg1_mask(ind_fit),0.95);
    bn0 = sqrt(area_max*dx*dy/pi);
    Cn0 = min(vg1_mask(ind_fit));
    X0n0 = interp1(1:size(X1,2),X1(1,:),indc(1));
    Y0n0 = interp1(1:size(Y1,1),Y1(:,1),indc(2));
    beta02D=[An0 bn0 Cn0 X0n0 Y0n0];
    model2D='gaussian2D';
    [beta2D,r2D,J2D,CovB2D,MSE] = nlinfit([X1(ind_fit) Y1(ind_fit)],vg1_mask(ind_fit),model2D,beta02D,options);
    An=beta2D(1)+beta2D(3); bn=beta2D(2);
    X0n=beta2D(4); Y0n=beta2D(5);
    Cn=beta2D(3);
    dn = sqrt((X0n-X0n0).^2+(Y0n-Y0n0).^2);
    if dn>Ri
        X0n = X0n0;
        Y0n = Y0n0;
        An = An0;
        bn = bn0;
        MSE = 1;
    end
    r2n0=(X1(ind_fit)-X0n0).^2 + (Y1(ind_fit)-Y0n0).^2;
    y_2D = vg1_mask(ind_fit);
    r2n=(X1(ind_fit)-X0n).^2 + (Y1(ind_fit)-Y0n).^2;
    nlin_y_2D=An.*exp(-r2n./(bn^2))+Cn;
    area_out(count) = area_max*dx*dy;
    x0_out(count) = X0n;
    y0_out(count) = Y0n;
    z0_out(count) = zlevel;
    b_out(count) = bn;
    I_out(count) = An*vg1_base;
    MSE_out(count) = MSE;
    
    % plot diagnostics
    if diagnosis == 1
        if mod(count-1,1)==0
            figure(3)
            subplot(1,2,1)
            pcolorjw(X1,Y1,10*log10(vg1));
            axis image
            %colormap('jet');
            ct = cbrewer('seq','Reds',9);
            colormap(ct)
            caxis([-90 -30]);
            xlim(xrange);
            ylim(yrange);
            hold on;
            plot(X0n,Y0n,'.g','markersize',20);
            plot(X0n0,Y0n0,'.b','markersize',20);
            hold off;
            xlabel('Easting of COVIS ( m )');
            ylabel('Northing of COVIS ( m )');
            h = colorbar;
            title(h,'vss (dB)');
            title(sprintf('height: %3.1f m',zlevel-z_vent));
            set(gca,'fontsize',12);
            
            subplot(1,2,2)
            pcolorjw(X1,Y1,10*log10(vg1_mask*vg1_base));
            axis image
            %colormap('jet');
            ct = cbrewer('seq','Reds',9);
            colormap(ct)
            caxis([-90 -30]);
            xlim(xrange);
            ylim(yrange);
            hold on;
            plot(X0n,Y0n,'.g','markersize',20);
            plot(X0n0,Y0n0,'.b','markersize',20);
            hold off;
            xlabel('Easting of COVIS ( m )');
            ylabel('Northing of COVIS ( m )');
            h = colorbar;
            title(h,'vss (dB)');
            title(sprintf('height: %3.1f m',zlevel-z_vent));
            set(gca,'fontsize',12);
                       
        end
        pause;
    end
    
    
    
    
end
center.area = area_out;
center.I = I_out;
center.b = b_out;
center.x = x0_out;
center.y = y0_out;
center.z = z0_out;
center.mse = MSE_out;
%end





