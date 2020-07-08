% This function calculates the plume radius, peak vertical velocity,
% background vertical velocity, and the coordinates of peak vertical
% velocity for each horizontal cross-section. 2-D Gaussian fit is performed
% on each cross-section to calculate these properties. The result is saved
% in the structural variable doppler_W.

function [doppler_W]=doppler_velocity(covis,Ig,zmin,zmax,xrange)

xg=covis.grid{1}.x;
yg=covis.grid{1}.y;
zg=covis.grid{1}.z;

ht_inc = covis.grid{1}.spacing.dz;
spacings = [0.5 0.5 0.5];
seta = 0.1;
[zs,i0,iN]=pick_heights(zg,zmin,zmax,ht_inc);



if isempty(xrange)
    xrange = [min(xg(:)),max(xg(:))];
end
yrange=[min(yg(:)),max(yg(:))];
Avsht = zeros(iN,3);
bvsht = zeros(iN,3);
Cvsht = zeros(iN,3);
for i=i0:iN
    % set z to process
    zlevel=zs(i);
    % calculate the vertical index
    z_min=min(min(min(zg)));
    zindex=floor((zlevel-z_min)./spacings(3))+1;
    xx=covis.grid{1}.bounds.xmin:spacings(1):covis.grid{2}.bounds.xmax;
    yy=covis.grid{1}.bounds.ymin:spacings(2):covis.grid{2}.bounds.ymax;
    ix=find(xx>=xrange(1)&xx<=xrange(2));
    iy=find(yy>=yrange(1)&yy<=yrange(2));
    % calculate the horizontal index
    [X1,Y1]=meshgrid(xx(ix),yy(iy));
    I1=Ig(iy,ix,zindex); %Vertical Velocity (m/s)
    
    
    Ibkg=median(I1(I1(:)>=0)); % calculate the background level (median)
    
    
    % Locate the maximum within each profile
    Imax=max(I1(:));
    [ii,jj]=find(I1==Imax);
    if Imax == 0||length(ii)~=1
        disp('meet a bad profile')
        Avsht(i,:)=[nan nan nan];
        bvsht(i,:)=[nan nan nan];
        Cvsht(i,:)=[nan nan nan];
        X0vsht(i,:)=[nan nan nan];
        Y0vsht(i,:)=[nan nan nan];
        bpred(i)=nan;
        R2vsht(i)=nan;
        normvsht(i)=nan;
        SSR(i)=nan; % sum_of_squares regression
        SSDev(i)=nan; % sum of squares residuals (deviation)
        SST(i)=nan; % sum of squares total
        Npts(i)=nan;
        continue;
    end
    x0=X1(ii,jj);
    y0=Y1(ii,jj);
    x0 = median(x0(:));
    y0 = median(y0(:));
    % extract the 1D profiles
    [pro,coor]=extract_profiles(x0,y0,X1,Y1,I1,spacings);
    i_x=find(pro{:,1}>=0);
    i_y=find(pro{:,2}>=0);
    if length(i_x)>=2&&length(i_y)>=2
        pro{1}(:)=interp1(coor{1}(i_x),pro{1}(i_x),coor{1}(:));
        pro{2}(:)=interp1(coor{2}(i_y),pro{2}(i_y),coor{2}(:));
    else
        Avsht(i,:)=[nan nan nan];
        bvsht(i,:)=[nan nan nan];
        Cvsht(i,:)=[nan nan nan];
        X0vsht(i,:)=[nan nan nan];
        Y0vsht(i,:)=[nan nan nan];
        bpred(i)=nan;
        R2vsht(i)=nan;
        normvsht(i)=nan;
        SSR(i)=nan; % sum_of_squares regression
        SSDev(i)=nan; % sum of squares residuals (deviation)
        SST(i)=nan; % sum of squares total
        Npts(i)=nan;
        continue;
    end
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
    An=beta2D(1); bn=beta2D(2)*1e6;
    X0n=beta2D(4)*1e6; Y0n=beta2D(5)*1e6;
    Cn=beta2D(3);
    r2n=(X1-X0n).^2 + (Y1-Y0n).^2;
    nlin_y_2D=An.*exp(-r2n./(bn^2))+Cn;
    
    % printing values
    %     fprintf('2D Gauss Fit for %s at height %d m\n',runname,zlevel)
    %     fprintf('Nonlinear Fit: C=(%4.2e)*exp(-r^2/(%5.2f)^2) + (%4.2e)\n',An,bn,Cn);
    %     fprintf('r^2=(X1-(%5.2f)).^2 + (Y1-(%5.2f)).^2;\n',X0n,Y0n)
    SS_res=sqrt(sum(sum(r2D.^2))/length(r2D(:))); % SSD
    %fprintf('Sum of squares (residuals) for nonlinear fit: %4.2e\n',SS_res);
    mean_nlin_y_2D = mean(nlin_y_2D(:));
    SS_reg=sum((nlin_y_2D(:)-mean_nlin_y_2D).^2);  % SSR
    %fprintf('Sum of squares (regression): %4.2e\n',SS_reg);
    SS_tot=sum((longC-mean_nlin_y_2D).^2);   % SST
    %fprintf('Sum of squares (total): %4.2e\n',SS_tot);
    R2=SS_reg./SS_tot;
    %fprintf('multiple correlation coeficient R^2=%5.2f  R=%5.2f\n',R2,sqrt(R2))
    [nrows,ncols]=size(nlin_y_2D);
    npts=nrows*ncols;
    normfit=norm(nlin_y_2D(:)-mean_nlin_y_2D);
    %fprintf('norm of residuals %5.2e\n',normfit)
    betaci=nlparci(beta2D,r2D,J2D);
    %fprintf('An: %5.2e < %5.2e < %5.2e', betaci(1,1),An,betaci(1,2))
    %fprintf('  vs. guessed %5.2e\n',Imax-Ibkg)
    %fprintf('bn: %5.2f < %5.2f < %5.2f', betaci(2,1)*1e6,bn,betaci(2,2)*1e6)
    %fprintf('  vs. predicted %5.2f\n',seta*zlevel)
    %fprintf('Cn: %5.2e < %5.2e < %5.2e', betaci(3,1),Cn,betaci(3,2))
    %fprintf('  vs. guessed %5.2e\n',Ibkg)
    %fprintf('X0n: %5.2f < %5.2f < %5.2f', betaci(4,1)*1e6,X0n,betaci(4,2)*1e6)
    %fprintf('  vs. guessed %5.2f\n',x0)
    %fprintf('Y0n: %5.2f < %5.2f < %5.2f', betaci(5,1)*1e6,Y0n,betaci(5,2)*1e6)
    %fprintf('  vs. guessed %5.2f\n',y0)
    
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



% save the results in a structural variable doppler_W
doppler_W.x=xx(ix);
doppler_W.y=yy(iy);
doppler_W.z=zs;
doppler_W.peak=Avsht;
doppler_W.radius=bvsht;
doppler_W.back=Cvsht;
doppler_W.x0=X0vsht;
doppler_W.y0=Y0vsht;
doppler_W.r = R2vsht;

