% This function calculates the plume radius, peak backscatter intensity,
% background intensity, and the coordinates of the peak backscatter
% intensity on each horizontal cross-section. 2-D Gaussian fit is performed
% on each cross-section to calculate those properties. The result is saved
% in the structure array doppler_I.

function [doppler_I]=doppler_intensity(covis,Ig,zmin,zmax,xrange)

if strcmp(covis.type,'doppler');
    xg=covis.grid{2}.x;
    yg=covis.grid{2}.y;
    zg=covis.grid{2}.z;
    dz=covis.grid{2}.spacing.dz;
    dx=covis.grid{2}.spacing.dx;
    dy=covis.grid{2}.spacing.dy;
    xx = covis.grid{2}.axis(1):dx:covis.grid{2}.axis(2);
    yy = covis.grid{2}.axis(3):dy:covis.grid{2}.axis(4);
    spacings = [covis.grid{2}.spacing.dx,covis.grid{2}.spacing.dy,covis.grid{2}.spacing.dz];
elseif strcmp(covis.type,'imaging');
    Ig(isnan(Ig))=0;
    xg=covis.grid.x;
    yg=covis.grid.y;
    zg=covis.grid.z;
    dz=covis.grid.spacing.dz;
    dx=covis.grid.spacing.dx;
    dy=covis.grid.spacing.dy;
    xx = covis.grid.axis(1):dx:covis.grid.axis(2);
    yy = covis.grid.axis(3):dy:covis.grid.axis(4);
    spacings = [covis.grid.spacing.dx,covis.grid.spacing.dy,covis.grid.spacing.dz];
end

seta = 0.1;
[zs,i0,iN]=pick_heights(zg,zmin,zmax,dz);


yrange=[min(yg(:)),max(yg(:))];
if isempty(xrange)
    xrange = [min(xg(:)),max(xg(:))];
end
zindex=zeros(1,iN);
for i=i0:iN
    % set z to process
    zlevel=zs(i);
    % calculate the vertical index
    zmin=min(min(min(zg)));
    zindex(i)=floor((zlevel-zmin)./dz)+1;
    ix=find(xx>=xrange(1)&xx<=xrange(2));
    iy=find(yy>=yrange(1)&yy<=yrange(2));
    % calculate the horizontal index
    [X1,Y1]=meshgrid(xx(ix),yy(iy));
    I1=Ig(iy,ix,zindex(i)); % back scatter intensity (m^-1)
    
    Ibkg=median(I1(:)); % calculate the background level (median)
    
    % Locate the maximum within each profile
    Imax=max(I1(:));
    [ii,jj]=find(I1==Imax);
    x0=X1(ii,jj);
    y0=Y1(ii,jj);
    x0 = median(x0(:));
    y0 = median(y0(:));
    % extract the 1D profiles
    [pro,coor]=extract_profiles(x0,y0,X1,Y1,I1,spacings);
    
    % 1D gauss fit
    model='gaussian';
    beta0y=[(Imax-Ibkg) (seta*zlevel)*10^-6 Ibkg y0*10^-6];
    [betay,ry,J] = nlinfit_xgy(coor{2},pro{2},model,beta0y);
    nlin_y_y=betay(1).*exp(-((coor{2}-(1e6*betay(4))).^2)./((1e6*betay(2)).^2))+betay(3);
    beta0x=[(Imax-Ibkg) (seta*zlevel)*1e-6 Ibkg x0*1e-6];
    [betax,rx,J] = nlinfit_xgy(coor{1},pro{1},model,beta0x);
    nlin_y_x=betax(1).*exp(-((coor{1}-(1e6*betax(4))).^2)./((1e6*betax(2)).^2))+betax(3);
    
    
    % 2D Gauss fit
    An0=mean([betax(1) betay(1)]);
    bn0=mean([betax(2) betay(2)])*1e6;
    X0n0=betax(4)*1e6;
    Y0n0=betay(4)*1e6;
    Cn0=mean([betax(3) betay(3)]);
    beta02D=[An0 bn0*1e-6 Cn0 X0n0*1e-6 Y0n0*1e-6];
    longX=X1(:); longY=Y1(:); longC=I1(:);
    model2D='gaussian2D';
    [beta2D,r2D,J2D] = nlinfit_xgy([longX longY],longC,model2D,beta02D);
    An=beta2D(1)+beta2D(3); bn=beta2D(2)*1e6;
    X0n=beta2D(4)*1e6; Y0n=beta2D(5)*1e6;
    Cn=beta2D(3);
    r2n=(X1-X0n).^2 + (Y1-Y0n).^2;
    nlin_y_2D=An.*exp(-r2n./(bn^2))+Cn;
    SS_res=sum(sum(r2D.^2)); % SSD
    %     fprintf('Sum of squares (residuals) for nonlinear fit: %4.2e\n',SS_res);
    mean_nlin_y_2D = mean(nlin_y_2D(:));
    SS_reg=sum((nlin_y_2D(:)-mean_nlin_y_2D).^2);  % SSR
    %     fprintf('Sum of squares (regression): %4.2e\n',SS_reg);
    SS_tot=sum((longC-mean_nlin_y_2D).^2);   % SST
    %     fprintf('Sum of squares (total): %4.2e\n',SS_tot);
    R2=SS_reg./SS_tot;
    %     fprintf('multiple correlation coeficient R^2=%5.2f  R=%5.2f\n',R2,sqrt(R2))
    [nrows,ncols]=size(nlin_y_2D);
    npts=nrows*ncols;
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
    
    % saving values for param. vs. height plots
    Avsht(i,:)=[betaci(1,1) An betaci(1,2)];
    bvsht(i,:)=[betaci(2,1)*1e6 bn betaci(2,2)*1e6];
    Cvsht(i,:)=[betaci(3,1) Cn betaci(3,2)];
    X0vsht(i,:)=[betaci(4,1)*1e6 X0n betaci(4,2)*1e6];
    Y0vsht(i,:)=[betaci(5,1)*1e6 Y0n betaci(5,2)*1e6];
    bpred(i)=seta*zlevel;
    R2vsht(i)=R2;
    normvsht(i)=normfit;
    SSR(i)=SS_reg; % sum_of_squares regression
    SSDev(i)=SS_res; % sum of squares residuals (deviation)
    SST(i)=SS_tot; % sum of squares total
    Npts(i)=npts;
end



% save the results in a structural variable doppler_I
doppler_I.z=zs;
doppler_I.x=xx(ix);
doppler_I.y=yy(iy);
doppler_I.peak=Avsht; % peak back scatter intensity (m^-1)
doppler_I.radius=bvsht; % radius of the plume (m)
doppler_I.back=Cvsht; % background intensity
doppler_I.x0=X0vsht; % maximum coordinate x-axis
doppler_I.y0=Y0vsht; % maximum coordinate y-axis
doppler_I.zindex=zindex;



