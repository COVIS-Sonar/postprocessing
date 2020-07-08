% This function calculates the plume radius, peak backscatter intensity,
% background intensity, and the coordinates of the peak backscatter
% intensity on each horizontal cross-section. 2-D Gaussian fit is performed
% on each cross-section to calculate those properties for the plume above
% Inferno.

function center = get_plume_centerline_properties_fun(grd,xs,ys,Ri,Ro,diagnosis)
% Input:
% grd: gridded data
% xs: x-coordinate of sourve vent (e.g., Inferno) (m)
% ys: y-coordinate of source vent (m)
% Ri: search radius for locating maximum in volume scattering
% Ro: search radius for finding nonlinear fit of a Gaussian function to a
% horizontal cross-section of volume scattering
% diagnosis: toggle for turning on and off the diagnostic plots
% Output:
% center: Matlab structure that includes estimated plume properties and
% diagnostic parameters. 

%load imaging data
% data_dir = 'F:\COVIS\Axial\COVIS_data\processed\Imaging\2018\08';
% %data_name = 'COVIS-20180819T063002-imaging2_sensor.mat';
% data_name = 'COVIS-20180819T070002-imaging1_sensor.mat';
% load(fullfile(data_dir,data_name));

% Example
% xs = -13.27;
% ys = -0.77;
% Ro = 8;
% Ri = 4;
% diagnosis = 1;


% set up grid parameters
xrange = [xs-Ro,xs+Ro];
yrange = [ys-Ro,ys+Ro];
zmin = 2.5;
zmax = 14;
xg=grd.x;
yg=grd.y;
zg=grd.z;
xx = squeeze(xg(1,:,1));
yy = squeeze(yg(:,1,1));
zz = squeeze(zg(1,1,:));
spacings = [grd.spacing.dx,grd.spacing.dy,grd.spacing.dz];




% main loop
fac = 10^-7;
seta = 0.076;
[~,i0] = min(abs(zz-zmin));
[~,iN] = min(abs(zz-zmax));
count = 0;
A_out = nan(iN-i0+1,3);
b_out = nan(size(A_out));
C_out = nan(size(A_out));
X0_out = nan(size(A_out));
Y0_out = nan(size(A_out));
Z0_out = nan(iN-i0+1,1);
bpred_out = nan(size(Z0_out));
R2_out = nan(size(Z0_out));
norm_out = nan(size(Z0_out));
SSR_out  = nan(size(Z0_out));
SSDev_out = nan(size(Z0_out));
SST_out = nan(size(Z0_out));
Npts_out  = nan(size(Z0_out));
Flag_out = nan(size(Z0_out));
Ig = grd.Id;
for iz=i0:iN
    count = count+1;
    zlevel=zz(iz);
    ix=find(xx>=xrange(1)&xx<=xrange(2));
    iy=find(yy>=yrange(1)&yy<=yrange(2));
    % calculate the horizontal index
    [X1,Y1]=meshgrid(xx(ix),yy(iy));
    I1=Ig(iy,ix,iz); % back scatter intensity (m^-1)
    %Ibkg=nanmedian(I1(:)); % calculate the background level (median)
    Ibkg = quantile(I1(:),0.1);
    
    % Locate the maximum within each profile
    d = sqrt((X1(:)-xs).^2+(Y1(:)-ys).^2);
    Imax=nanmax(I1(d<Ri));
    if isnan(Imax)
        continue
    end
    [ii,jj]=find(I1==Imax);
    x0=X1(ii,jj);
    y0=Y1(ii,jj);
    x0 = median(x0(:));
    y0 = median(y0(:));
    % extract the 1D profiles
    [pro,coor]=extract_profiles(x0,y0,X1,Y1,I1,spacings);
    
    % 1D gauss fit
    model='gaussian';
    beta0y=[(Imax-Ibkg) (seta*zlevel)*fac Ibkg y0*fac];
    [betay,ry,J,flagx] = nlinfit_xgy(coor{2},pro{2},model,beta0y,fac);
    nlin_y_y=betay(1).*exp(-((coor{2}-(betay(4)/fac)).^2)./((betay(2)/fac).^2))+betay(3);
    beta0x=[(Imax-Ibkg) (seta*zlevel)*fac Ibkg x0*fac];
    [betax,rx,J,flagy] = nlinfit_xgy(coor{1},pro{1},model,beta0x,fac);
    nlin_y_x=betax(1).*exp(-((coor{1}-(betax(4)/fac)).^2)./((betax(2)/fac).^2))+betax(3);
    
    
    % 2D Gauss fit
    An0=mean([betax(1) betay(1)]);
    bn0=mean([betax(2) betay(2)])/fac;
    X0n0=betax(4)/fac;
    Y0n0=betay(4)/fac;
    %Cn0=mean([betax(3) betay(3)]);
    Cn0 = quantile(I1(~isnan(I1)),0.05);
    beta02D=[An0 bn0*fac Cn0 X0n0*fac Y0n0*fac];
    longX=X1(:); longY=Y1(:); longC=I1(:);
    model2D='gaussian2D';
    [beta2D,r2D,J2D,flagxy] = nlinfit_xgy([longX longY],longC,model2D,beta02D,fac);
    An=beta2D(1)+beta2D(3); bn=beta2D(2)/fac;
    X0n=beta2D(4)/fac; Y0n=beta2D(5)/fac;
    Cn=beta2D(3);
    r2n=(X1-X0n).^2 + (Y1-Y0n).^2;
    nlin_y_2D=An.*exp(-r2n./(bn^2))+Cn;
    SS_res=sum(sum(r2D.^2)); % SSD
    %     fprintf('Sum of squares (residuals) for nonlinear fit: %4.2e\n',SS_res);
    mean_nlin_y_2D = mean(nlin_y_2D(:));
    SS_reg=sum((nlin_y_2D(:)-mean_nlin_y_2D).^2);  % SSR
    %     fprintf('Sum of squares (regression): %4.2e\n',SS_reg);
    SS_tot=sum((longC(~isnan(longC))-mean_nlin_y_2D).^2);   % SST
    %     fprintf('Sum of squares (total): %4.2e\n',SS_tot);
    R2=SS_reg./SS_tot;
    %     fprintf('multiple correlation coeficient R^2=%5.2f  R=%5.2f\n',R2,sqrt(R2))
    npts = length(longC(~isnan(longC)));
    normfit=norm(nlin_y_2D(:)-mean_nlin_y_2D);
    %     fprintf('norm of residuals %5.2e\n',normfit)
    betaci=nlparci(beta2D,r2D,J2D);
    %     fprintf('An: %5.2e < %5.2e < %5.2e', betaci(1,1),An,betaci(1,2))
    % 	    fprintf('  vs. guessed %5.2e\n',Imax-Ibkg)
    %     fprintf('bn: %5.2f < %5.2f < %5.2f', betaci(2,1)*1e6,bn,betaci(2,2)*1e6)
    % 	    fprintf('  vs. predicted %5.2f\n',seta*zlevel)
    %     fprintf('Cn: %5.2e < %5.2e < %5.2e', betaci(3,1),Cn,betaci(3,2))
    % 	    fprintf('  vs. guessed %5.2e\n',Ibkg)
    %     fprintf('X0n: %5.2f < %5.2f < %5.2f', betaci(4,1)*1e6,X0n,betaci(4,2)*1e6)
    % 	    fprintf('  vs. guessed %5.2f\n',x0)
    %     fprintf('Y0n: %5.2f < %5.2f < %5.2f', betaci(5,1)*1e6,Y0n,betaci(5,2)*1e6)
    %     	fprintf('  vs. guessed %5.2f\n',y0)
    
    % plot diagnostics
    if diagnosis == 1
        if mod(count,4)==0
            figure
            subplot(1,2,1)
            plot(sqrt(r2n(:)),longC(:),'.');
            hold on;
            plot(sqrt(r2n(:)),nlin_y_2D(:),'r');
            hold off;
            xlabel('Distance ( m )')
            ylabel('dB');
            %xlim([0 10]);
            %ylim([-90 -30]);
            title(sprintf('height: %f',zlevel));
            subplot(1,2,2)
            pcolorjw(xx,yy,10*log10(Ig(:,:,iz)));
            caxis([-90 -30]);
            axis image
            colormap('jet');
            xlim(xrange);
            ylim(yrange);
            hold on;
            plot(X0n,Y0n,'.k','markersize',10);
            hold off;
            xlabel('Easting of COVIS ( m )');
            ylabel('Northing of COVIS ( m )');
            h = colorbar;
            title(h,'dB');
            title(sprintf('height: %f',zlevel));
        end
    end
    if flagx+flagy+flagxy>0
        flag = 1;
    else
        flag = 0;
    end
    
    % save values for param. vs. height plots
    A_out(count,:)=[betaci(1,1) An betaci(1,2)];
    b_out(count,:)=[betaci(2,1)*fac bn betaci(2,2)*fac];
    C_out(count,:)=[betaci(3,1) Cn betaci(3,2)];
    X0_out(count,:)=[betaci(4,1)*fac X0n betaci(4,2)*fac];
    Y0_out(count,:)=[betaci(5,1)*fac Y0n betaci(5,2)*fac];
    Z0_out(count) = zlevel;
    bpred_out(count)=seta*zlevel;
    R2_out(count)=R2;
    norm_out(count)=normfit;
    SSR_out(count)=SS_reg; % sum_of_squares regression
    SSDev_out(count)=SS_res; % sum of squares residuals (deviation)
    SST_out(count)=SS_tot; % sum of squares total
    Npts_out(count)=npts;
    Flag_out(count) = flag;
end
center.I = A_out;
center.be = b_out;
center.I_bkg = C_out;
center.x = X0_out;
center.y = Y0_out;
center.z = Z0_out;
center.R2 = R2_out;
center.norm = norm_out;
center.SSR = SSR_out;
center.SSDev = SSDev_out;
center.SST = SST_out;
center.Npts = Npts_out;
center.Flag = Flag_out;
end





