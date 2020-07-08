% this function is used to filter the output axial velocity data using the
% window function and connected component algorithm.
function [v_filt,v_std_filt,area,vr_back,vz_back] = doppler_filter(v,v_std,vr,Ig,doppler_I,index,method)
% calculate the threshold used to filter the velocity output bay assuming a Gaussian distribution for the VSS
% profile at each horizontal corss-section. The threshold is set to be the
% VSS value at a distance R times the radius away from the centerline. 


R = 2;
R_b = 4;
windthresh=doppler_I.peak(index,2)*exp(-R^2)+doppler_I.back(index,2);
windthresh_back = doppler_I.peak(index,2)*exp(-7)+doppler_I.back(index,2);
% calculate the background velocity
win= exp(-(Ig/(windthresh_back)).^4);
vr_back = vr.*win;
vz_back = v.*win;
switch lower(method)
    case 'window'
    x = -40:0.5:-5;
    y = -25:0.5:5;
    [xx,yy]=meshgrid(x,y);
    c_x = doppler_I.x0(index,2);
    c_y = doppler_I.y0(index,2);
    b_I=doppler_I.radius(index,2); % radius of the intensity profile (m)
    rr = sqrt((xx-c_x).^2+(yy-c_y).^2);
    i_r = find(rr>R_b*b_I); % the radius of the circle is 3 times the Gaussian radius
    v(i_r) = 0;
    v_std(i_r) = 0;
    win= 1 - exp(-(Ig/(windthresh)).^4);
    v_filt=v.*win;
    v_std_filt=v_std.*win;

    case 'ccomp'

    s = Ig>windthresh;
    ww = zeros(size(v));
    ww(s)=1;
    ww(~s)=0;
    cc = bwconncomp(ww);
    %object_num = cc.NumObjects;
    %if object_num >=1
        object_prop = regionprops(cc);
        object_area = [object_prop.Area];
        [area_max,i_area]=max(object_area);
        ind = cc.PixelIdxList;
        ww_d = ww;
        ww_d(ind{i_area})=0;
        ww = ww-ww_d;
    %end
    v_filt = v.*ww;
    v_std_filt = v_std.*ww;
    I_filt = Ig.*ww;
    area = area_max*0.5*0.5;
    if isempty(area)
        area = nan;
    end
end