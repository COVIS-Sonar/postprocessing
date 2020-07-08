function [pro,coor]=extract_profiles(x0,y0,X1,Y1,I1,spacings)
% modified by Guangyu Xu July 2012
pro = cell(1,2);
coor = cell(1,2);
xmin=min(min(X1));
xindex=floor((x0-xmin)./spacings(1))+1; 
ymin=min(min(Y1));
yindex=floor((y0-ymin)./spacings(2))+1; 
pro{1}=I1(yindex,:)';
coor{1}=X1(yindex,:)';
pro{2}=I1(:,xindex);
coor{2}=Y1(:,xindex);
%[ylen,xlen]=size(I1);


%prolen=min([xlen,ylen]);

% if xindex<(prolen./2),
%     xstart=1;
%     xend=prolen;
%     %fprintf('pro 1 (%d,%d:%d)\n',yindex,xstart,xend);
% else
%     xstart=xlen-prolen+1;
%     xend=xlen;
%     %fprintf('pro 1 (%d,%d:%d)\n',yindex,xstart,xend);
% end
% if yindex<(prolen./2),
%     ystart=1;
%     yend=prolen;
%     %fprintf('pro 2 (%d:%d,%d)\n',ystart,yend,xindex);
% else
%     ystart=ylen-prolen+1;
%     yend=ylen;
%     %fprintf('pro 2 (%d:%d,%d)\n',ystart,yend,xindex);
% end


% pro(:,1)=I1(yindex,xstart:xend)';
% coor(:,1)=X1(yindex,xstart:xend)';
% pro(:,2)=I1(ystart:yend,xindex);
% coor(:,2)=Y1(ystart:yend,xindex);
