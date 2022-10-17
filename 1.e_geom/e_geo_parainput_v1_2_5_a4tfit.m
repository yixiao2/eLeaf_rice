% eLeaf: 3D model of rice leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 1.2.6

function e_geo_parainput_v1_2_5_a4tfit(CFG_PARA_COM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0 for CKIR64; 1 for HCIR64
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%select_com=input('selection (e.g [1 1 0 0 0 0 0 0 0]):');
select_com=CFG_PARA_COM;
if size(select_com,2)~=9
    disp('input selection error');
    return;
end

GLB_digits=10;% digits to keep in COMSOL and MATLAB. Also used when passing values to 'parameters' in COMSOL
% Added to keep robustness of ray tracing.

if select_com(1)==0
    MST_thickatBU=31.40326721e-6;%1st
    MST_thickatvein=64.42342204e-6;%2nd
    LEAF_thickatvein=75.51998995e-6;%3rd
    LEAF_thickatBU=70.49649136e-6;%4th
    VEIN_width=186.4871879e-6;%7th
else
    MST_thickatBU=27.90968e-6;%1st
    MST_thickatvein=67.41458e-6;%2nd
    LEAF_thickatvein=77.98927e-6;%3rd
    LEAF_thickatBU=71.30234e-6;%4th
    VEIN_width=200.4024e-6;%7th
end
MST_thickatBU=str2num(num2str(MST_thickatBU,GLB_digits));% later use same GLB_digits when passing to 'parameters' in COMSOL
MST_thickatvein=str2num(num2str(MST_thickatvein,GLB_digits));
LEAF_thickatvein=str2num(num2str(LEAF_thickatvein,GLB_digits));
LEAF_thickatBU=str2num(num2str(LEAF_thickatBU,GLB_digits));
VEIN_width=str2num(num2str(VEIN_width,GLB_digits));

if select_com(2)==0
    BS_min_width=9.6540668e-6;%5th
    BS_area=870.2022e-12;%6th
    BS_plastid_volume=0.1044446;% unit %
else
    BS_min_width=9.0126916e-6;%5th
    BS_area=657.5816e-12;%6th
    BS_plastid_volume=0.1176085;
end

%EPL_thick=EPU_thick=(Leaf Thickness-Mesophyll Thickness)/2
EPL_thick=(LEAF_thickatvein-MST_thickatvein)/2;
EPU_thick=EPL_thick;
%BU_thick=Leaf Thickness at bulliform-mesophyll thickness at bulliform-EPL_thick
BU_thick=LEAF_thickatBU-MST_thickatBU-EPL_thick;
BS_thick=BS_min_width;
VEIN_r=(BS_area/pi-BS_min_width^2)/(2*BS_min_width);
%VEIN_r=double(vpa(VEIN_r,4));
if VEIN_r<0 || (VEIN_r+BS_min_width)*2>MST_thickatvein
    display('Please check leaf level anatomical parameters');
    %    break;
end
BSE_width=VEIN_r*2;
EPU_width=(VEIN_width-BSE_width)/2/2;
BUT_width=(VEIN_width-BSE_width)/2;
VEIN_l=LEAF_thickatvein/2;
%BS_chlo_thick=2e-6; x^2-tmp_B*x+tmp_C=0
BS_dthick=1e-6;
tmp_C=BS_plastid_volume*((VEIN_r+BS_thick)^2-VEIN_r^2);
tmp_B=2*(VEIN_r+BS_thick-BS_dthick);
BS_chlo_thick=(tmp_B-sqrt(tmp_B^2-4*tmp_C))/2;

EPL_thick=str2num(num2str(EPL_thick,GLB_digits));
EPU_thick=str2num(num2str(EPU_thick,GLB_digits));
BU_thick=str2num(num2str(BU_thick,GLB_digits));
BS_thick=str2num(num2str(BS_thick,GLB_digits));
VEIN_r=str2num(num2str(VEIN_r,GLB_digits));
BSE_width=str2num(num2str(BSE_width,GLB_digits));
EPU_width=str2num(num2str(EPU_width,GLB_digits));
BUT_width=str2num(num2str(BUT_width,GLB_digits));
VEIN_l=str2num(num2str(VEIN_l,GLB_digits));
BS_chlo_thick=str2num(num2str(BS_chlo_thick,GLB_digits));

if select_com(3)==0
    IAS_3D_input=0.0802;% 8th; relative to whole leaf in 3D
else
    %% replace 2.porosity
    IAS_3D_input=0.0617;% 8th; mesophyll porosity
    %% end replace
end

if select_com(8)==0
    chl_con=0.000543;% 9th; chl content mg/mm2
    chl_con_ratio_BvM=1;% ratio
else
    %% replace 3.chlorophyll content
    chl_con=0.000551;% 9th; chl content mg/mm2
    chl_con_ratio_BvM=1;
    %% end replace
end

if select_com(4)==0
    cell_length=18.98022162;%10th
    cell_height=13.25452238;%11th
    cell_volume=1886.56;%12th
    mit_r=0.05;
    mit_l=0.05;
    vac_lx=0.05;
    vac_ly=0.05;
    epsln=0.2;
    epsln4=0;
    flatness=3;
    lmbd=0.2;rho=0.12;
    %distance for rho: 1/sqrt(3)*rho*(scale_x+scale_y)/2
else
    %% replace 4.cell shape
    cell_length=19.77811;%10th
    cell_height=13.35439;%11th
    cell_volume=1908.5;%12th
    mit_r=0.05;
    mit_l=0.05;
    vac_lx=0.05;
    vac_ly=0.05;
    epsln=0.2;
    epsln4=0;
    flatness=3;
    lmbd=0.2;rho=0.08;
    %distance for rho: 1/sqrt(3)*rho*(scale_x+scale_y)/2
    %% end replace 4
end
cell_length=str2num(num2str(cell_length,GLB_digits));
cell_height=str2num(num2str(cell_height,GLB_digits));
%cell_thick in e_geo_MS3D.m

if select_com(7)==0
    SmS=0.316251875;% 13th; CW touching another cell / total external CW length
else
    %% replace 5.Sm
    SmS=0.420127845;% 13th; CW touching another cell / total external CW length
    %% end replace 5
end

if select_com(6)==0
    %for e_physics
    %cellwallthick=0.16;% 14th
    cellwallthick=0.217565;
else
    %% replace 6.wallthick
    %cellwallthick=0.17414763;%14th
    cellwallthick=0.229958;
    %% end replace 6
end

if select_com(5)==0
    plastid_volume=0.6971006;%15th
else
    %% replace 7.plastid vol
    plastid_volume=0.6438284;%15th
    %% end replace
end

if select_com(9)==0
    Resp=2.92;
    yiill=0.63;
    s_yin=0.63;
    Jm=224.63;
    theta=0.98;%constant
    Vcmax=114.39;
    Sco=3375;%constant
    Kc=239;
    Ko=266000;%Km=239*(1+210000/266000)=428
    
    fc_CA=1.26;%fold change of CA
    fc_wall_perm=2.27;%fold change of wall permeability
    
    %bst_ft_vm=1.4;
    %bst_ft_jm=1.2;
    %bst_ft_rd=-0.2;
    %bst_ft_s=1.0;
    %bst_ft_yii=0.9;
    %Vcmax=Vcmax*bst_ft_vm;
    %Jm=Jm*bst_ft_jm;
    %Resp=Resp*bst_ft_rd;
    %s_yin=s_yin*bst_ft_s;
    %yiill=yiill*bst_ft_yii;
else
    %% replace 8.enzyme
    Resp=3.34;
    yiill=0.61;
    s_yin=0.48;
    Jm=268.79;
    theta=0.98;%constant
    Vcmax=121.98;
    Sco=3375;%constant
    Kc=239;
    Ko=266000;%Km=239*(1+210000/266000)=428
    %% end replace 8
    
    fc_CA=1.26;%fold change of CA
    fc_wall_perm=2.27;%fold change of wall permeability
    
    %bst_ft_vm=3.3;
    %bst_ft_jm=1.7;
    %bst_ft_rd=0.9;
    %bst_ft_s=1.5;
    %bst_ft_yii=0.95;
    %Vcmax=Vcmax*bst_ft_vm;
    %Jm=Jm*bst_ft_jm;
    %Resp=Resp*bst_ft_rd;
    %s_yin=s_yin*bst_ft_s;
    %yiill=yiill*bst_ft_yii;
end

MS_dthick=1e-6;
BS_dthick=1e-6;

BS_mito_l=2e-6;
BS_mito_r=0.5e-6;

Nlobe=8;
dw=0.05;
%dse=0.5;
dvaz=0.1;

save parainput.mat
display('Input parameter loaded.')
end
