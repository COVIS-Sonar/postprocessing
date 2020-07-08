function [zs,i0,iN]=pick_heights(zg,min_start_ht,max_end_ht,ht_inc)
%   function [zs,i0,iN]=pick_heights(zg,min_start_ht,max_end_ht,ht_inc)
%
%
%
minht=min(zg(:)); 
%fprintf('minht is %f \n',minht)
maxht=max(zg(:));
%fprintf('maxht is %f \n',maxht)
if minht<min_start_ht,
    minht=min_start_ht;
else
   diffht=minht-min_start_ht;
   steps_up=ceil(diffht./ht_inc);
   minht=min_start_ht+steps_up*ht_inc;
end
if maxht>max_end_ht,
    maxht=max_end_ht;
else
   diffht=max_end_ht-maxht;
   steps_down=ceil(diffht./ht_inc);
   maxht=max_end_ht-steps_down*ht_inc;
end
zs=minht:ht_inc:maxht;
i0=1;
iN=max(size(zs));
fprintf('Height coverage is %d m to %d m.\n',minht,maxht);
fprintf('Indices are %d to %d.\n',i0,iN);

