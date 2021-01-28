% eLeaf: 3D model of rice leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 1.2.5

%MST_thickatBU=31.40326721e-6;%1st
%MST_thickatvein=64.42342204e-6;%2nd
%LEAF_thickatvein=75.51998995e-6;%3rd
    %EPL_thick=EPU_thick=(Leaf Thickness-Mesophyll Thickness)/2
    EPL_thick=(LEAF_thickatvein-MST_thickatvein)/2;
    %EPL_thick=double(vpa(EPL_thick,4));% 4 digits
    EPU_thick=EPL_thick;
%LEAF_thickatBU=70.49649136e-6;%4th
    %BU_thick=Leaf Thickness at bulliform-mesophyll thickness at bulliform-EPL_thick
    BU_thick=LEAF_thickatBU-MST_thickatBU-EPL_thick;
%BS_min_width=6.183979e-6;%5th
    BS_thick=BS_min_width;
%BS_area=304.6302e-12;%6th
    VEIN_r=(BS_area/pi-BS_min_width^2)/(2*BS_min_width);
    %VEIN_r=double(vpa(VEIN_r,4));
    
    if VEIN_r<0 || (VEIN_r+BS_min_width)*2>MST_thickatvein
        display('Please check leaf level anatomical parameters');
        %    break;
    end
    
    BSE_width=VEIN_r*2;
%VEIN_width=186.4871879e-6;%7th
    EPU_width=(VEIN_width-BSE_width)/2/2;
    BUT_width=(VEIN_width-BSE_width)/2;
    VEIN_l=LEAF_thickatvein/2;
