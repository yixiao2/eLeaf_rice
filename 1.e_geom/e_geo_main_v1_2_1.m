% eLeaf: 3D model of rice leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 1.2.6

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update 2017/12/27
% - mph geometry pass the test of free tetrahedron mesh
% before export object for ray tracing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% WARNING: this specific version fix dse in line22
loop_idx=-1;
loop_idy=-1;
status=0;

clearvars -except ms_distribution N_frame loop_frame tmpnum
tic;tic;
loop_idx=-1;loop_idy=-1;status=0;
load parainput.mat
load tmp_MS3D.mat lmbd rho cell_thick dse leaf_z0

import com.comsol.model.*
import com.comsol.model.util.*

% try

model = ModelUtil.create('Model_3Dleaf');

%BSE_width=20e-6;
model.param.set('BSE_width', [num2str(BSE_width),'[m]']);
%EPL_thick=10e-6;
model.param.set('EPL_thick', [num2str(EPL_thick),'[m]'], 'thickness of lower epidermis');
%MST_thickatvein=65e-6;
model.param.set('MST_thickatvein', [num2str(MST_thickatvein),'[m]'], 'mesophyll thickness at vein');
%EPU_width=20e-6;
model.param.set('EPU_width', [num2str(EPU_width),'[m]']);
%BUT_width=120e-6;
model.param.set('BUT_width', [num2str(BUT_width),'[m]']);
model.param.set('MST_thickatBU', [num2str(MST_thickatBU),'[m]']);
model.param.set('VEIN_r', [num2str(VEIN_r),'[m]']);
model.param.set('BS_thick', [num2str(BS_thick),'[m]']);
model.param.set('VEIN_l', [num2str(VEIN_l),'[m]']);
model.param.set('EPU_thick', [num2str(EPU_thick),'[m]']);
model.param.set('BU_thick', [num2str(BU_thick),'[m]']);
model.param.set('MS_dthick', [num2str(MS_dthick),'[m]']);
model.param.set('BS_dthick', [num2str(BS_dthick),'[m]']);
model.param.set('BS_chlo_thick', [num2str(BS_chlo_thick),'[m]']);
%LEAF_z=cell_thick/2*1.1*1e-6;
LEAF_z=leaf_z0;
model.param.set('LEAF_z', [num2str(LEAF_z),'[m]']);
BS_chlo_z=cell_thick/2*0.9*1e-6;
model.param.set('BS_chlo_z', [num2str(BS_chlo_z),'[m]']);
model.param.set('BS_mito_l', [num2str(BS_mito_l),'[m]']);
model.param.set('BS_mito_r', [num2str(BS_mito_r),'[m]']);

model.modelNode.create('mod1');

model.geom.create('geom1', 3);
%model.geom('geom1').geomRep('comsol');
%% DEFAULT relative repair tolerance
%model.geom('geom1').repairTol(1.0E-8);
model.geom('geom1').repairTolType('auto');
model.geom('geom1').feature.create('wp1', 'WorkPlane');
model.geom('geom1').feature.create('ext1', 'Extrude');
model.geom('geom1').feature.create('ext2', 'Extrude');
model.geom('geom1').feature.create('sph1', 'Sphere');
model.geom('geom1').feature.create('sph2', 'Sphere');
model.geom('geom1').feature.create('sph3', 'Sphere');
model.geom('geom1').feature.create('sph4', 'Sphere');
model.geom('geom1').feature.create('sph5', 'Sphere');
model.geom('geom1').feature.create('sph6', 'Sphere');
model.geom('geom1').feature.create('uni1', 'Union');
model.geom('geom1').feature.create('int1', 'Intersection');
model.geom('geom1').feature.create('dif1', 'Difference');
model.geom('geom1').feature.create('del1', 'Delete');
model.geom('geom1').feature('wp1').geom.feature.create('b1', 'BezierPolygon');
model.geom('geom1').feature('wp1').geom.feature.create('pol1', 'Polygon');
model.geom('geom1').feature('wp1').geom.feature.create('c1', 'Circle');
model.geom('geom1').feature('wp1').geom.feature.create('b2', 'BezierPolygon');
model.geom('geom1').feature('wp1').geom.feature.create('c2', 'Circle');
model.geom('geom1').feature('wp1').geom.feature.create('b3', 'BezierPolygon');
model.geom('geom1').feature('wp1').geom.feature.create('c3', 'Circle');
model.geom('geom1').feature('wp1').geom.feature.create('c4', 'Circle');
model.geom('geom1').feature('wp1').geom.feature.create('c5', 'Circle');
model.geom('geom1').feature('wp1').geom.feature.create('dif1', 'Difference');
model.geom('geom1').feature('wp1').geom.feature.create('dif2', 'Difference');
model.geom('geom1').feature('wp1').geom.feature.create('dif3', 'Difference');
model.geom('geom1').feature('wp1').geom.feature.create('dif4', 'Difference');
model.geom('geom1').feature('wp1').geom.feature.create('dif5', 'Difference');
model.geom('geom1').feature('wp1').geom.feature.create('dif6', 'Difference');
model.geom('geom1').feature('wp1').geom.feature.create('del1', 'Delete');
model.geom('geom1').feature('wp1').geom.feature('b1').set('p', {'BSE_width/2' 'BSE_width/2' 'BSE_width/2+EPU_width' 'BSE_width/2+EPU_width+BUT_width/6' 'BSE_width/2+EPU_width+BUT_width/6*2' 'BSE_width/2+EPU_width+BUT_width/6*3' 'BSE_width/2+EPU_width+BUT_width/6*3'; 'EPL_thick' 'EPL_thick+MST_thickatvein' 'EPL_thick+MST_thickatvein' 'EPL_thick+MST_thickatvein' 'EPL_thick+MST_thickatBU' 'EPL_thick+MST_thickatBU' 'EPL_thick'});
model.geom('geom1').feature('wp1').geom.feature('b1').set('degree', {'1' '1' '3' '1'});
model.geom('geom1').feature('wp1').geom.feature('b1').set('w', {'1' '1' '1' '1' '1' '1' '1' '1' '1' '1'});
model.geom('geom1').feature('wp1').geom.feature('pol1').set('source', 'table');
model.geom('geom1').feature('wp1').geom.feature('pol1').set('table', {'0' 'EPL_thick'; '0' 'EPL_thick+MST_thickatvein'; 'BSE_width/2' 'EPL_thick+MST_thickatvein'; 'BSE_width/2' 'EPL_thick'});
model.geom('geom1').feature('wp1').geom.feature('c1').set('angle', '180');
model.geom('geom1').feature('wp1').geom.feature('c1').set('rot', '-90');
model.geom('geom1').feature('wp1').geom.feature('c1').set('r', 'VEIN_r+BS_thick');
model.geom('geom1').feature('wp1').geom.feature('c1').set('pos', {'0' 'VEIN_l'});
model.geom('geom1').feature('wp1').geom.feature('b2').set('p', {'0' '0' 'BSE_width/2+EPU_width' 'BSE_width/2+EPU_width+BUT_width/6' 'BSE_width/2+EPU_width+BUT_width/6*2' 'BSE_width/2+EPU_width+BUT_width/6*3' 'BSE_width/2+EPU_width+BUT_width/6*3'; '0' 'EPL_thick+MST_thickatvein+EPU_thick' 'EPL_thick+MST_thickatvein+EPU_thick' 'EPL_thick+MST_thickatvein+EPU_thick' 'EPL_thick+MST_thickatBU+BU_thick' 'EPL_thick+MST_thickatBU+BU_thick' '0'});
model.geom('geom1').feature('wp1').geom.feature('b2').set('degree', {'1' '1' '3' '1'});
model.geom('geom1').feature('wp1').geom.feature('b2').set('w', {'1' '1' '1' '1' '1' '1' '1' '1' '1' '1'});
model.geom('geom1').feature('wp1').geom.feature('c2').name('VEIN');
model.geom('geom1').feature('wp1').geom.feature('c2').set('angle', '180');
model.geom('geom1').feature('wp1').geom.feature('c2').set('rot', '-90');
model.geom('geom1').feature('wp1').geom.feature('c2').set('r', 'VEIN_r');
model.geom('geom1').feature('wp1').geom.feature('c2').set('pos', {'0' 'VEIN_l'});
model.geom('geom1').feature('wp1').geom.feature('b3').set('p', {'BSE_width/2+MS_dthick' 'BSE_width/2+MS_dthick' 'BSE_width/2+EPU_width' 'BSE_width/2+EPU_width+BUT_width/6' 'BSE_width/2+EPU_width+BUT_width/6*2' 'BSE_width/2+EPU_width+BUT_width/6*3-+MS_dthick' 'BSE_width/2+EPU_width+BUT_width/6*3-MS_dthick'; 'EPL_thick+MS_dthick' 'EPL_thick+MST_thickatvein-MS_dthick' 'EPL_thick+MST_thickatvein-MS_dthick' 'EPL_thick+MST_thickatvein-MS_dthick' 'EPL_thick+MST_thickatBU-MS_dthick' 'EPL_thick+MST_thickatBU-MS_dthick' 'EPL_thick+MS_dthick'});
model.geom('geom1').feature('wp1').geom.feature('b3').set('w', {'1' '1' '1' '1' '1' '1' '1' '1' '1' '1'});
model.geom('geom1').feature('wp1').geom.feature('b3').set('degree', {'1' '1' '3' '1'});
model.geom('geom1').feature('wp1').geom.feature('c3').set('r', 'VEIN_r+BS_thick+MS_dthick');
model.geom('geom1').feature('wp1').geom.feature('c3').set('angle', '180');
model.geom('geom1').feature('wp1').geom.feature('c3').set('pos', {'0' 'VEIN_l'});
model.geom('geom1').feature('wp1').geom.feature('c3').set('rot', '-90');
model.geom('geom1').feature('wp1').geom.feature('c4').set('r', 'VEIN_r+BS_thick-BS_dthick');
model.geom('geom1').feature('wp1').geom.feature('c4').set('angle', '180');
model.geom('geom1').feature('wp1').geom.feature('c4').set('pos', {'0' 'VEIN_l'});
model.geom('geom1').feature('wp1').geom.feature('c4').set('rot', '-90');
model.geom('geom1').feature('wp1').geom.feature('c5').set('r', 'VEIN_r+BS_thick-BS_dthick-BS_chlo_thick');
model.geom('geom1').feature('wp1').geom.feature('c5').set('angle', '180');
model.geom('geom1').feature('wp1').geom.feature('c5').set('pos', {'0' 'VEIN_l'});
model.geom('geom1').feature('wp1').geom.feature('c5').set('rot', '-90');
model.geom('geom1').feature('wp1').geom.feature('dif1').name('BS');
model.geom('geom1').feature('wp1').geom.feature('dif1').set('keep', true);
model.geom('geom1').feature('wp1').geom.feature('dif1').set('intbnd', false);
model.geom('geom1').feature('wp1').geom.feature('dif1').set('edge', 'all');
model.geom('geom1').feature('wp1').geom.feature('dif1').selection('input2').set({'c2'});
model.geom('geom1').feature('wp1').geom.feature('dif1').selection('input').set({'c1'});
model.geom('geom1').feature('wp1').geom.feature('dif2').name('BSE');
model.geom('geom1').feature('wp1').geom.feature('dif2').set('keep', true);
model.geom('geom1').feature('wp1').geom.feature('dif2').set('intbnd', false);
model.geom('geom1').feature('wp1').geom.feature('dif2').set('edge', 'all');
model.geom('geom1').feature('wp1').geom.feature('dif2').selection('input2').set({'c1'});
model.geom('geom1').feature('wp1').geom.feature('dif2').selection('input').set({'pol1'});
model.geom('geom1').feature('wp1').geom.feature('dif3').name('EP');
model.geom('geom1').feature('wp1').geom.feature('dif3').set('keep', true);
model.geom('geom1').feature('wp1').geom.feature('dif3').set('intbnd', false);
model.geom('geom1').feature('wp1').geom.feature('dif3').set('edge', 'all');
model.geom('geom1').feature('wp1').geom.feature('dif3').selection('input2').set({'b1' 'pol1'});
model.geom('geom1').feature('wp1').geom.feature('dif3').selection('input').set({'b2'});
model.geom('geom1').feature('wp1').geom.feature('dif4').name('MS');
model.geom('geom1').feature('wp1').geom.feature('dif4').set('keep', true);
model.geom('geom1').feature('wp1').geom.feature('dif4').set('intbnd', false);
model.geom('geom1').feature('wp1').geom.feature('dif4').set('edge', 'all');
model.geom('geom1').feature('wp1').geom.feature('dif4').selection('input2').set({'c1'});
model.geom('geom1').feature('wp1').geom.feature('dif4').selection('input').set({'b1'});
model.geom('geom1').feature('wp1').geom.feature('dif5').name('MS_i');
model.geom('geom1').feature('wp1').geom.feature('dif5').set('intbnd', false);
model.geom('geom1').feature('wp1').geom.feature('dif5').set('edge', 'all');
model.geom('geom1').feature('wp1').geom.feature('dif5').selection('input2').set({'c3'});
model.geom('geom1').feature('wp1').geom.feature('dif5').selection('input').set({'b3'});
model.geom('geom1').feature('wp1').geom.feature('dif6').name('BS_chlo');
model.geom('geom1').feature('wp1').geom.feature('dif6').set('intbnd', false);
model.geom('geom1').feature('wp1').geom.feature('dif6').set('edge', 'all');
model.geom('geom1').feature('wp1').geom.feature('dif6').selection('input2').set({'c5'});
model.geom('geom1').feature('wp1').geom.feature('dif6').selection('input').set({'c4'});
model.geom('geom1').feature('wp1').geom.feature('del1').selection('input').init;
model.geom('geom1').feature('wp1').geom.feature('del1').selection('input').set({'b1' 'b2' 'c1' 'pol1'});
model.geom('geom1').feature('ext1').setIndex('distance', 'LEAF_z', 0);
model.geom('geom1').feature('ext1').set('face', 'all');
model.geom('geom1').feature('ext1').set('crossfaces', false);
model.geom('geom1').feature('ext1').selection('input').set({'wp1.c2' 'wp1.dif1' 'wp1.dif2' 'wp1.dif3' 'wp1.dif4' 'wp1.dif5'});
%geom('geom1').feature('ext2').set('createselection', true);%save tag
count_list_cresel=0;
count_list_cresel=count_list_cresel+1;
list_cresel(count_list_cresel)=cellstr('ext2');
prep_CreSel(count_list_cresel,1)=count_list_cresel;
prep_CreSel(count_list_cresel,2)=9;%9 means BS chlo
model.geom('geom1').feature('ext2').setIndex('distance', 'BS_chlo_z', 0);
model.geom('geom1').feature('ext2').selection('input').set({'wp1.dif6'});
tmp_str=['bs_chl'];
model.geom('geom1').feature('ext2').name(tmp_str);
model.geom('geom1').feature('sph1').set('r', 'BS_mito_r');
model.geom('geom1').feature('sph1').set('pos', {'(VEIN_r+BS_thick-BS_dthick-BS_chlo_thick-BS_mito_l)*sin(15/180*pi)' '(VEIN_r+BS_thick-BS_dthick-BS_chlo_thick-BS_mito_l)*cos(15/180*pi)+VEIN_l' '0'});
model.geom('geom1').feature('sph2').set('r', 'BS_mito_r');
model.geom('geom1').feature('sph2').set('pos', {'(VEIN_r+BS_thick-BS_dthick-BS_chlo_thick-BS_mito_l)*sin((15+30)/180*pi)' '(VEIN_r+BS_thick-BS_dthick-BS_chlo_thick-BS_mito_l)*cos((15+30)/180*pi)+VEIN_l' '0'});
model.geom('geom1').feature('sph3').set('r', 'BS_mito_r');
model.geom('geom1').feature('sph3').set('pos', {'(VEIN_r+BS_thick-BS_dthick-BS_chlo_thick-BS_mito_l)*sin((15+30*2)/180*pi)' '(VEIN_r+BS_thick-BS_dthick-BS_chlo_thick-BS_mito_l)*cos((15+30*2)/180*pi)+VEIN_l' '0'});
model.geom('geom1').feature('sph4').set('r', 'BS_mito_r');
model.geom('geom1').feature('sph4').set('pos', {'(VEIN_r+BS_thick-BS_dthick-BS_chlo_thick-BS_mito_l)*sin((15+30*3)/180*pi)' '(VEIN_r+BS_thick-BS_dthick-BS_chlo_thick-BS_mito_l)*cos((15+30*3)/180*pi)+VEIN_l' '0'});
model.geom('geom1').feature('sph5').set('r', 'BS_mito_r');
model.geom('geom1').feature('sph5').set('pos', {'(VEIN_r+BS_thick-BS_dthick-BS_chlo_thick-BS_mito_l)*sin((15+30*4)/180*pi)' '(VEIN_r+BS_thick-BS_dthick-BS_chlo_thick-BS_mito_l)*cos((15+30*4)/180*pi)+VEIN_l' '0'});
model.geom('geom1').feature('sph6').set('r', 'BS_mito_r');
model.geom('geom1').feature('sph6').set('pos', {'(VEIN_r+BS_thick-BS_dthick-BS_chlo_thick-BS_mito_l)*sin((15+30*5)/180*pi)' '(VEIN_r+BS_thick-BS_dthick-BS_chlo_thick-BS_mito_l)*cos((15+30*5)/180*pi)+VEIN_l' '0'});
model.geom('geom1').feature('uni1').selection('input').set({'sph1' 'sph2' 'sph3' 'sph4' 'sph5' 'sph6'});
model.geom('geom1').feature('dif1').name('BS_mito');
%model.geom('geom1').feature('int1').set('createselection', true);%save tag
count_list_cresel=count_list_cresel+1;
list_cresel(count_list_cresel)=cellstr('int1');
prep_CreSel(count_list_cresel,1)=count_list_cresel;
prep_CreSel(count_list_cresel,2)=10;
model.geom('geom1').feature('int1').set('keep', true);
model.geom('geom1').feature('int1').set('intbnd', false);
model.geom('geom1').feature('int1').set('face', 'all');
model.geom('geom1').feature('int1').selection('input').set({'uni1' 'ext1(2)'});
tmp_str=['bs_mits'];
model.geom('geom1').feature('int1').name(tmp_str);
model.geom('geom1').feature('dif1').name('bs_cyto');
%model.geom('geom1').feature('dif1').set('createselection', true);%save tag
count_list_cresel=count_list_cresel+1;
list_cresel(count_list_cresel)=cellstr('dif1');
prep_CreSel(count_list_cresel,1)=count_list_cresel;
prep_CreSel(count_list_cresel,2)=8;
model.geom('geom1').feature('dif1').set('keep', true);
model.geom('geom1').feature('dif1').set('intbnd', false);
model.geom('geom1').feature('dif1').set('face', 'all');
model.geom('geom1').feature('dif1').selection('input2').set({'ext2' 'uni1'});
model.geom('geom1').feature('dif1').selection('input').set({'ext1(2)'});
model.geom('geom1').run('dif1');
% leaf section volume
model.geom('geom1').measure.selection.init(3);
tmp_Ndomains=model.geom('geom1').obj('ext1(1)').getNDomains;
eval(['model.geom(''geom1'').measure.selection.set(''ext1(1)'',',['[',num2str(1:tmp_Ndomains),']'],');'])
vol_vein=model.geom('geom1').measure().getVolume();
model.geom('geom1').measure.selection.init(3);
tmp_Ndomains=model.geom('geom1').obj('ext1(2)').getNDomains;
eval(['model.geom(''geom1'').measure.selection.set(''ext1(2)'',',['[',num2str(1:tmp_Ndomains),']'],');'])
vol_bs=model.geom('geom1').measure().getVolume();
model.geom('geom1').measure.selection.init(3);
tmp_Ndomains=model.geom('geom1').obj('ext1(3)').getNDomains;
eval(['model.geom(''geom1'').measure.selection.set(''ext1(3)'',',['[',num2str(1:tmp_Ndomains),']'],');'])
vol_bse=model.geom('geom1').measure().getVolume();
model.geom('geom1').measure.selection.init(3);
tmp_Ndomains=model.geom('geom1').obj('ext1(4)').getNDomains;
eval(['model.geom(''geom1'').measure.selection.set(''ext1(4)'',',['[',num2str(1:tmp_Ndomains),']'],');'])
vol_ep=model.geom('geom1').measure().getVolume();
model.geom('geom1').measure.selection.init(3);
tmp_Ndomains=model.geom('geom1').obj('ext1(5)').getNDomains;
eval(['model.geom(''geom1'').measure.selection.set(''ext1(5)'',',['[',num2str(1:tmp_Ndomains),']'],');'])
vol_msregion=model.geom('geom1').measure().getVolume();

model.geom('geom1').feature('del1').selection('input').init;
model.geom('geom1').feature('del1').selection('input').set({'uni1' 'ext1(2)'});
model.geom('geom1').run('del1');


%% Mesophyll shape
%Nlobe=8;% number of lobes; optional 4+2n (n>=0)
model.param.set('Nlobe',num2str(Nlobe));
Tlobe=Nlobe/2;% type of lobe
model.param.set('Tlobe',num2str(Tlobe));
%dw=0.05;
model.param.set('dw', num2str(dw));
%dse=0.5;
model.param.set('dse', num2str(dse));
model.param.set('dvaz', num2str(dvaz));
%cell_length=14.7;%unit um
%cell_height=10.6;%unit um
MSC_length=cell_length*1e-6;
MSC_height=cell_height*1e-6;
%cell_volume=1783;%unit um^3
%calculate cell_thick in z axis
%cell_thick=round(cell_volume*3/4/pi/cell_length*2/cell_height*2*2,3,'significant');
model.param.set('cell_length',[num2str(cell_length),'[um]']);
model.param.set('cell_height',[num2str(cell_height),'[um]']);
model.param.set('cell_thick',[num2str(cell_thick),'[um]']);
%mit_r=0.07;
%mit_l=0.05;
model.param.set('mit_r',num2str(mit_r));
model.param.set('mit_l',num2str(mit_l));
model.param.set('vac_lx',num2str(vac_lx));
model.param.set('vac_ly',num2str(vac_ly));
model.param.set('lmbd',num2str(lmbd));
model.param.set('rho',num2str(rho));
%epsln=0.2;
if Nlobe==4
    model.param.set('epsln','0');
else
    model.param.set('epsln',num2str(epsln));
end
model.param.set('epsln4',num2str(epsln4));%special epsln for cell with 4 lobes
model.param.set('scale_x','(cell_length/2)/(1/2*Tlobe)');
model.param.set('scale_y','(cell_height/2)/(1/(2*sqrt(3))*(Tlobe+1))');
tmp_Npts=14+4*(Tlobe-2)+2*Nlobe;

%% generate basic mesophyll cell with different lobes
count_wp=2;
count_ext=3;
count_elp=1;
count_int=2;
count_dif=2;
count_mov=1;
if loop_idx==-1&&loop_idy==-1
    % random the initial point
    %     Origin_x=BSE_width/2-MSC_length/2*rand;
    %     Origin_y=EPL_thick-MSC_height/2*rand;
    load tmp_IAS2D.mat Origin_x Origin_y
    %delete('save_tmp0.mat')
else
    % adjust Origin_x and Origin_y according to loop_idx and loop_idy
    load save_tmp0.mat Origin_x Origin_y
    max_loop_num=11;
    seed_dx=0.1*MSC_length/max_loop_num*(loop_idx-(max_loop_num-1)/2);
    seed_dy=0.1*MSC_height/max_loop_num*(loop_idy-(max_loop_num-1)/2);
    Origin_x=Origin_x+seed_dx;
    Origin_y=Origin_y+seed_dy;
end

for loop_lobe=0:Tlobe-2
    tmp_Nlobe=Nlobe-2*loop_lobe;
    tmp_Tlobe=tmp_Nlobe/2;
    if tmp_Nlobe==4
        tmp_epsln=epsln4;
    else
        tmp_epsln=epsln;
    end
    [tmp_wallString,tmp_chloString,tmp_chliString,tmp_IDXxse,tmp_IDXyse]=msoutline_v1_2(tmp_Nlobe,dw,dse,tmp_epsln);
    tmptag_wp=['wp',num2str(count_wp)];count_wp=count_wp+1;
    model.geom('geom1').feature.create(tmptag_wp, 'WorkPlane');
    model.geom('geom1').feature(tmptag_wp).geom.feature.create('pol1', 'Polygon');
    model.geom('geom1').feature(tmptag_wp).geom.feature.create('pol2', 'Polygon');
    model.geom('geom1').feature(tmptag_wp).geom.feature.create('pol3', 'Polygon');
    model.geom('geom1').feature(tmptag_wp).set('quickz', '-cell_thick');
    model.geom('geom1').feature(tmptag_wp).geom.feature('pol1').set('source', 'table');
    eval(['model.geom(''geom1'').feature(tmptag_wp).geom.feature(''pol1'').set(''table'', {',tmp_wallString,'});']);
    model.geom('geom1').feature(tmptag_wp).geom.feature('pol2').set('source', 'table');
    eval(['model.geom(''geom1'').feature(tmptag_wp).geom.feature(''pol2'').set(''table'', {',tmp_chloString,'});']);
    model.geom('geom1').feature(tmptag_wp).geom.feature('pol3').set('source', 'table');
    eval(['model.geom(''geom1'').feature(tmptag_wp).geom.feature(''pol3'').set(''table'', {',tmp_chliString,'});']);
    tmptag_ext=['ext',num2str(count_ext)];count_ext=count_ext+1;
    model.geom('geom1').feature.create(tmptag_ext, 'Extrude');
    model.geom('geom1').feature(tmptag_ext).setIndex('distance', 'cell_thick*2', 0);
    model.geom('geom1').feature(tmptag_ext).selection('input').set(tmptag_wp);
    rab=1;
    ea=rab*sqrt((tmp_Tlobe/2)^2/rab^2+25/4/3);
    eb=sqrt((tmp_Tlobe/2)^2/rab^2+25/4/3);
    tmptag_elp=['elp',num2str(count_elp)];count_elp=count_elp+1;
    model.geom('geom1').feature.create(tmptag_elp, 'Ellipsoid');
    eval(['model.geom(''geom1'').feature(tmptag_elp).set(''semiaxes'', {''scale_x*',num2str(ea),''' ''scale_y*',num2str(eb),''' ''cell_thick/2''});'])
    tmptag_int=['int',num2str(count_int)];count_int=count_int+1;
    tag_sca2=tmptag_int;
    model.geom('geom1').feature.create(tmptag_int, 'Intersection');
    model.geom('geom1').feature(tmptag_int).set('face', 'all');
    model.geom('geom1').feature(tmptag_int).set('intbnd', false);
    model.geom('geom1').feature(tmptag_int).selection('input').set([cellstr(tmptag_elp),cellstr([tmptag_ext,'(1)'])]);
    tmptag_elp=['elp',num2str(count_elp)];count_elp=count_elp+1;
    model.geom('geom1').feature.create(tmptag_elp, 'Ellipsoid');
    eval(['model.geom(''geom1'').feature(tmptag_elp).set(''semiaxes'', {''scale_x*(',num2str(ea),'-dw)'' ''scale_y*(',num2str(eb),'-dw)'' ''cell_thick/2-dw*(scale_x+scale_y)/2''});'])
    tmptag_int=['int',num2str(count_int)];count_int=count_int+1;
    tag_sca3=tmptag_int;
    model.geom('geom1').feature.create(tmptag_int, 'Intersection');
    model.geom('geom1').feature(tmptag_int).set('face', 'all');
    model.geom('geom1').feature(tmptag_int).set('intbnd', false);
    model.geom('geom1').feature(tmptag_int).selection('input').set([cellstr(tmptag_elp),cellstr([tmptag_ext,'(2)'])]);
    tmptag_elp=['elp',num2str(count_elp)];count_elp=count_elp+1;
    model.geom('geom1').feature.create(tmptag_elp, 'Ellipsoid');
    if tmp_Nlobe==4%% chl_i
        eval(['model.geom(''geom1'').feature(tmptag_elp).set(''semiaxes'', {''scale_x*(',num2str(ea),'-dse*(1-dw+epsln4/2))'' ''scale_y*(',num2str(eb),'-dse*(1-dw+epsln4/2))'' ''cell_thick/2-dse*(1-dw+epsln4/2)*(scale_x+scale_y)/2''});'])
    else
        eval(['model.geom(''geom1'').feature(tmptag_elp).set(''semiaxes'', {''scale_x*(',num2str(ea),'-dse*(1-dw+epsln/2))'' ''scale_y*(',num2str(eb),'-dse*(1-dw+epsln/2))'' ''cell_thick/2-dse*(1-dw+epsln/2)*(scale_x+scale_y)/2''});'])
    end
    tmptag_int=['int',num2str(count_int)];count_int=count_int+1;
    tag_sca4=tmptag_int;
    model.geom('geom1').feature.create(tmptag_int, 'Intersection');
    model.geom('geom1').feature(tmptag_int).set('face', 'all');
    model.geom('geom1').feature(tmptag_int).set('intbnd', false);
    model.geom('geom1').feature(tmptag_int).selection('input').set([cellstr(tmptag_elp),cellstr([tmptag_ext,'(3)'])]);
    tmptag_elp=['elp',num2str(count_elp)];count_elp=count_elp+1;
    model.geom('geom1').feature.create(tmptag_elp, 'Ellipsoid');%% vacuole
    if tmp_Nlobe==4
        eval(['model.geom(''geom1'').feature(tmptag_elp).set(''semiaxes'', {''min(-scale_x*(',tmp_IDXxse{1},'+(mit_r+mit_l)/sqrt(3)*2+mit_r+vac_lx),-scale_x*(',tmp_IDXxse{3},'))'' ''min(scale_y*(',tmp_IDXyse{6},'-(mit_r+mit_l)*2-mit_r-vac_ly),scale_y*(',tmp_IDXyse{3},'))'' ''cell_thick/2-(dse*(1-dw+epsln/2)+dvaz)*(scale_x+scale_y)/2''});'])%x: pt6; y:pt8
    else
        eval(['model.geom(''geom1'').feature(tmptag_elp).set(''semiaxes'', {''-scale_x*(((',tmp_IDXxse{1},'+(mit_r+mit_l)/sqrt(3)*2)+1/2*(-',num2str(tmp_Tlobe),'+3))/2+mit_r+vac_lx)'' ''scale_y*(sqrt(3)/2*(1/2*(-',num2str(tmp_Tlobe),'+3)-(',tmp_IDXxse{1},'+(mit_r+mit_l)/sqrt(3)*2))-mit_r-vac_ly)'' ''cell_thick/2-(dse*(1-dw+epsln/2)+dvaz)*(scale_x+scale_y)/2''});'])%x: pt6; y:pt8
    end
    tag_sca1=tmptag_elp;
    %     if tmp_Nlobe==4
    %         model.geom('geom1').feature(tmptag_elp).active(false);
    %     end
    for tmp_loop=1:tmp_Nlobe %% mitochondria
        tmptag_elp=['elp',num2str(count_elp)];count_elp=count_elp+1;
        model.geom('geom1').feature.create(tmptag_elp, 'Ellipsoid');
        if tmp_loop==1
            tag_mits=cellstr(tmptag_elp);
        else
            tag_mits=[tag_mits,cellstr(tmptag_elp)];
        end
        model.geom('geom1').feature(tmptag_elp).set('semiaxes', {'scale_x*mit_r' 'scale_y*mit_r' 'cell_thick*mit_r'});
        if tmp_loop==1
            eval(['model.geom(''geom1'').feature(tmptag_elp).set(''pos'', {''scale_x*(',tmp_IDXxse{1},'+(mit_r+mit_l)/sqrt(3)*2)'' ''0'' ''0''});'])
        elseif tmp_loop==2
            if tmp_Nlobe==4
                eval(['model.geom(''geom1'').feature(tmptag_elp).set(''pos'', {''scale_x*(',num2str((-tmp_Tlobe/2-1+tmp_loop)),')'' ''scale_y*((',tmp_IDXyse{6},')-(mit_r+mit_l)*2)'' ''0''});'])
            else
                eval(['model.geom(''geom1'').feature(tmptag_elp).set(''pos'', {''scale_x*((',tmp_IDXxse{1},'+(mit_r+mit_l)/sqrt(3)*2)+1/2*(-',num2str(tmp_Tlobe),'+3))/2'' ''scale_y*sqrt(3)/2*(1/2*(-',num2str(tmp_Tlobe),'+3)-(',tmp_IDXxse{1},'+(mit_r+mit_l)/sqrt(3)*2))'' ''0''});'])
            end
        elseif tmp_loop>2&&tmp_loop<tmp_Nlobe/2
            eval(['model.geom(''geom1'').feature(tmptag_elp).set(''pos'', {''scale_x*(',num2str((-tmp_Tlobe/2-1+tmp_loop)),')'' ''scale_y*sqrt(3)/2*(1/2*(-',num2str(tmp_Tlobe),'+3)-(',tmp_IDXxse{1},'+(mit_r+mit_l)/sqrt(3)*2))'' ''0''});'])
        elseif tmp_loop==tmp_Nlobe/2
            eval(['model.geom(''geom1'').feature(tmptag_elp).set(''pos'', {''-scale_x*((',tmp_IDXxse{1},'+(mit_r+mit_l)/sqrt(3)*2)+1/2*(-',num2str(tmp_Tlobe),'+3))/2'' ''scale_y*sqrt(3)/2*(1/2*(-',num2str(tmp_Tlobe),'+3)-(',tmp_IDXxse{1},'+(mit_r+mit_l)/sqrt(3)*2))'' ''0''});'])
        elseif tmp_loop==tmp_Nlobe/2+1
            eval(['model.geom(''geom1'').feature(tmptag_elp).set(''pos'', {''-scale_x*(',tmp_IDXxse{1},'+(mit_r+mit_l)/sqrt(3)*2)'' ''0'' ''0''});'])
        elseif tmp_loop==tmp_Nlobe/2+2
            if tmp_Nlobe==4
                eval(['model.geom(''geom1'').feature(tmptag_elp).set(''pos'', {''scale_x*(',num2str((tmp_Tlobe/2+1-(tmp_loop-tmp_Tlobe))),')'' ''-scale_y*((',tmp_IDXyse{6},')-(mit_r+mit_l)*2)'' ''0''});'])
            else
                eval(['model.geom(''geom1'').feature(tmptag_elp).set(''pos'', {''-scale_x*((',tmp_IDXxse{1},'+(mit_r+mit_l)/sqrt(3)*2)+1/2*(-',num2str(tmp_Tlobe),'+3))/2'' ''-scale_y*sqrt(3)/2*(1/2*(-',num2str(tmp_Tlobe),'+3)-(',tmp_IDXxse{1},'+(mit_r+mit_l)/sqrt(3)*2))'' ''0''});'])
            end
        elseif tmp_loop>tmp_Nlobe/2+2&&tmp_loop<tmp_Nlobe
            eval(['model.geom(''geom1'').feature(tmptag_elp).set(''pos'', {''scale_x*(',num2str((tmp_Tlobe/2+1-(tmp_loop-tmp_Tlobe))),')'' ''-scale_y*sqrt(3)/2*(1/2*(-',num2str(tmp_Tlobe),'+3)-(',tmp_IDXxse{1},'+(mit_r+mit_l)/sqrt(3)*2))'' ''0''});'])
        else %tmp_loop==Nlobe
            eval(['model.geom(''geom1'').feature(tmptag_elp).set(''pos'', {''scale_x*((',tmp_IDXxse{1},'+(mit_r+mit_l)/sqrt(3)*2)+1/2*(-',num2str(tmp_Tlobe),'+3))/2'' ''-scale_y*sqrt(3)/2*(1/2*(-',num2str(tmp_Tlobe),'+3)-(',tmp_IDXxse{1},'+(mit_r+mit_l)/sqrt(3)*2))'' ''0''});'])
        end
    end
    tmptag_dif=['dif',num2str(count_dif)];count_dif=count_dif+1;
    model.geom('geom1').feature.create(tmptag_dif, 'Difference');
    model.geom('geom1').feature(tmptag_dif).set('face', 'all');
    model.geom('geom1').feature(tmptag_dif).set('intbnd', false);
    model.geom('geom1').feature(tmptag_dif).set('keep', false);
    model.geom('geom1').feature(tmptag_dif).selection('input').set(tag_sca3);
    model.geom('geom1').feature(tmptag_dif).selection('input2').set(tag_sca4);
    tmptag_mov=['mov',num2str(count_mov)];count_mov=count_mov+1;tag_types(loop_lobe+1,1)=cellstr(tmptag_mov);
    model.geom('geom1').feature.create(tmptag_mov, 'Move');
    model.geom('geom1').feature(tmptag_mov).selection('input').set(tag_sca1);
    model.geom('geom1').feature(tmptag_mov).set('displx',[num2str(Origin_x),'-0.5*scale_x*(',num2str(Tlobe-tmp_Tlobe),')']);%7.5e-6%3.3e-6[FAIL]
    model.geom('geom1').feature(tmptag_mov).set('disply',num2str(Origin_y));%6.9e-6%6.6e-6[FAIL]
    tmptag_mov2=['mov',num2str(count_mov)];count_mov=count_mov+1;tag_types(loop_lobe+1,2)=cellstr(tmptag_mov2);
    model.geom('geom1').feature.duplicate(tmptag_mov2, tmptag_mov);
    model.geom('geom1').feature(tmptag_mov2).selection('input').set(tmptag_dif);
    tmptag_mov2=['mov',num2str(count_mov)];count_mov=count_mov+1;tag_types(loop_lobe+1,3)=cellstr(tmptag_mov2);
    model.geom('geom1').feature.duplicate(tmptag_mov2, tmptag_mov);
    model.geom('geom1').feature(tmptag_mov2).selection('input').set(tag_mits);
    tmptag_mov2=['mov',num2str(count_mov)];count_mov=count_mov+1;tag_types(loop_lobe+1,4)=cellstr(tmptag_mov2);
    model.geom('geom1').feature.duplicate(tmptag_mov2, tmptag_mov);
    model.geom('geom1').feature(tmptag_mov2).selection('input').set(tag_sca2);
    model.geom('geom1').run(tmptag_mov2);
end

%% columns and rows
N_ms_cols=floor((EPU_width+BUT_width/2)/MSC_length)+2;
N_ms_rows=floor(MST_thickatvein/MSC_height)+2;

model.param.set('col_dx', '(Nlobe/2-1/2)*scale_x');
model.param.set('col_dy', '3/2/sqrt(3)*scale_y');
model.param.set('row_dx', '1/2*scale_x');
model.param.set('row_dy', '9/2/sqrt(3)*scale_y');

count_pseudoMSC=0;
model.geom('geom1').measure.selection.init(3);
model.geom('geom1').measure.selection.set('mov1', [1]);
vol_vac0=model.geom('geom1').measure().getVolume();
model.geom('geom1').measure.selection.init(3);
model.geom('geom1').measure.selection.set('mov3(1)', [1]);
vol_mit0=model.geom('geom1').measure().getVolume();
count_copy=1;
% count_int=5;
% count_dif=3;
count_del=2;
count_uni=2;
count_MSC=1;
load tmp_IAS2D.mat ms_distribution
% load 0000.mat ms_distribution
for loop_i=1:N_ms_rows
    for loop_j=0:N_ms_cols
        if ms_distribution(loop_i,loop_j+1,1)~=0&&ms_distribution(loop_i,loop_j+1,2)~=Tlobe
            if ms_distribution(loop_i,loop_j+1,1)==2
                ms_type=ms_distribution(loop_i,loop_j+1,2)+1;
                tmp_Nlobe=Nlobe-ms_distribution(loop_i,loop_j+1,2)*2;
                tmp_dx=ms_distribution(loop_i,loop_j+1,3);
                %%%%%%%%%%%%debug%%%%%%%%%%%%%
                %ms_type=1;
                %tmp_Nlobe=Nlobe;
                %tmp_dx=0;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else
                ms_type=1;
                tmp_Nlobe=Nlobe;
                tmp_dx=0;
            end
            %%%%%%%%%%%%%%%%%%%debug%%%%%%%%%%%%%%%%
            %         if count_MSC>35%Maximum 35 MSCs
            %             break;
            %         end
            %%%%%%%%%%%%%%%%%%%debug%%%%%%%%%%%%%%%%%
            if loop_j==0&&loop_i==1
                count_pseudoMSC=count_pseudoMSC+1;
            else
                %copy vacuole&mito&chlo&whole cell object
                %%copy mov1(vacuole)
                tmptag_vac0=['copy',num2str(count_copy)];%0 means the original; 1 means original intersect with the ms_i region
                model.geom('geom1').feature.create(tmptag_vac0, 'Copy');
                count_copy=count_copy+1;
                eval(['model.geom(''geom1'').feature(tmptag_vac0).set(''displx'',''col_dx*',num2str(loop_j),'-floor(',num2str(loop_j),'/3)*row_dx+(',num2str(loop_i),'-1)*row_dx+scale_x*1*(',num2str(tmp_dx),')'');']);
                eval(['model.geom(''geom1'').feature(tmptag_vac0).set(''disply'',''col_dy*',num2str(loop_j),'-floor(',num2str(loop_j),'/3)*row_dy+(',num2str(loop_i),'-1)*row_dy'');']);
                %model.geom('geom1').feature(tmptag_vac0).selection('input').set({'mov1'});
                model.geom('geom1').feature(tmptag_vac0).selection('input').set(tag_types{ms_type,1});
                model.geom('geom1').run(tmptag_vac0);
                tmptag_vac1=['int',num2str(count_int)];
                model.geom('geom1').feature.create(tmptag_vac1, 'Intersection');
                count_int=count_int+1;
                model.geom('geom1').feature(tmptag_vac1).set('keep', 'on');
                model.geom('geom1').feature(tmptag_vac1).set('intbnd', 'off');
                eval(['model.geom(''geom1'').feature(tmptag_vac1).selection(''input'').set({''',tmptag_vac0,''' ''ext1(6)''});']);
                %model.geom('geom1').feature(tmptag_vac1).set('createselection', true);%create selection before run
                model.geom('geom1').run(tmptag_vac1);
                tmp_getNDomains=model.geom('geom1').obj(tmptag_vac1).getNDomains;
                if tmp_getNDomains==0
                    %                 model.geom('geom1').feature(tmptag_vac0).active(false);
                    %                 model.geom('geom1').feature(tmptag_vac1).active(false);
                    model.geom('geom1').feature.remove(tmptag_vac0);
                    model.geom('geom1').feature.remove(tmptag_vac1);
                else
                    %%%intersect 10% vacuole
                    model.geom('geom1').measure.selection.init(3);
                    model.geom('geom1').measure.selection.set(tmptag_vac1, [1]);
                    tmp_vol_vac=model.geom('geom1').measure().getVolume();
                    if tmp_vol_vac/(vol_vac0/2)<0.01||tmp_vol_vac>vol_vac0%/2 because this threshold should be relative to half of the mov1(vacuole) %0.01 is 0.1 in v1
                        %                     model.geom('geom1').feature(tmptag_vac0).active(false);
                        %                     model.geom('geom1').feature(tmptag_vac1).active(false);
                        model.geom('geom1').feature.remove(tmptag_vac0);
                        model.geom('geom1').feature.remove(tmptag_vac1);
                    else
                        %%copy mov2(chlo)
                        tmptag_chl0=['copy',num2str(count_copy)];
                        model.geom('geom1').feature.create(tmptag_chl0, 'Copy');
                        count_copy=count_copy+1;
                        eval(['model.geom(''geom1'').feature(tmptag_chl0).set(''displx'',''col_dx*',num2str(loop_j),'-floor(',num2str(loop_j),'/3)*row_dx+(',num2str(loop_i),'-1)*row_dx+scale_x*1*(',num2str(tmp_dx),')'');']);
                        eval(['model.geom(''geom1'').feature(tmptag_chl0).set(''disply'',''col_dy*',num2str(loop_j),'-floor(',num2str(loop_j),'/3)*row_dy+(',num2str(loop_i),'-1)*row_dy'');']);
                        model.geom('geom1').feature(tmptag_chl0).selection('input').set(tag_types{ms_type,2});
                        tmptag_chl1=['int',num2str(count_int)];
                        model.geom('geom1').feature.create(tmptag_chl1, 'Intersection');
                        count_int=count_int+1;
                        model.geom('geom1').feature(tmptag_chl1).set('keep', 'on');
                        model.geom('geom1').feature(tmptag_chl1).set('intbnd', 'off');
                        eval(['model.geom(''geom1'').feature(tmptag_chl1).selection(''input'').set({''',tmptag_chl0,''' ''ext1(6)''});']);
                        %model.geom('geom1').feature(tmptag_chl1).set('createselection', true);%create selection before run
                        model.geom('geom1').run(tmptag_chl1);
                        %##operation int error
                        if model.geom('geom1').obj(tmptag_chl1).getNDomains==0
                            %                         model.geom('geom1').feature(tmptag_vac0).active(false);
                            %                         model.geom('geom1').feature(tmptag_vac1).active(false);
                            model.geom('geom1').feature.remove(tmptag_vac0);
                            model.geom('geom1').feature.remove(tmptag_vac1);
                            %                         model.geom('geom1').feature(tmptag_chl0).active(false);
                            %                         model.geom('geom1').feature(tmptag_chl1).active(false);
                            model.geom('geom1').feature.remove(tmptag_chl0);
                            model.geom('geom1').feature.remove(tmptag_chl1);
                        else
                            model.geom('geom1').measure.selection.init(3);
                            model.geom('geom1').measure.selection.set(tmptag_chl1, [1]);
                            tmp_vol1=model.geom('geom1').measure().getVolume();
                            model.geom('geom1').measure.selection.init(3);
                            model.geom('geom1').measure.selection.set(tmptag_chl0, [1]);
                            tmp_vol0=model.geom('geom1').measure().getVolume();
                            if tmp_vol1>tmp_vol0
                                %                             model.geom('geom1').feature(tmptag_vac0).active(false);
                                %                             model.geom('geom1').feature(tmptag_vac1).active(false);
                                model.geom('geom1').feature.remove(tmptag_vac0);
                                model.geom('geom1').feature.remove(tmptag_vac1);
                                %                             model.geom('geom1').feature(tmptag_chl0).active(false);
                                %                             model.geom('geom1').feature(tmptag_chl1).active(false);
                                model.geom('geom1').feature.remove(tmptag_chl0);
                                model.geom('geom1').feature.remove(tmptag_chl1);
                            else
                                %%copy mov3(mitos)
                                count_mit=0;
                                clear tmptag_mitall tmptag_mitall_del
                                for loop_k=1:tmp_Nlobe% remember the number of mitos remained
                                    tmptag_mit0=['copy',num2str(count_copy)];
                                    model.geom('geom1').feature.create(tmptag_mit0, 'Copy');
                                    count_copy=count_copy+1;
                                    eval(['model.geom(''geom1'').feature(tmptag_mit0).set(''displx'',''col_dx*',num2str(loop_j),'-floor(',num2str(loop_j),'/3)*row_dx+(',num2str(loop_i),'-1)*row_dx+scale_x*1*(',num2str(tmp_dx),')'');']);
                                    eval(['model.geom(''geom1'').feature(tmptag_mit0).set(''disply'',''col_dy*',num2str(loop_j),'-floor(',num2str(loop_j),'/3)*row_dy+(',num2str(loop_i),'-1)*row_dy'');']);
                                    eval(['model.geom(''geom1'').feature(tmptag_mit0).selection(''input'').set({''',char(tag_types{ms_type,3}),'(',num2str(loop_k),')''});']);
                                    tmptag_mit1=['int',num2str(count_int)];
                                    model.geom('geom1').feature.create(tmptag_mit1, 'Intersection');
                                    count_int=count_int+1;
                                    model.geom('geom1').feature(tmptag_mit1).set('keep', 'on');
                                    model.geom('geom1').feature(tmptag_mit1).set('intbnd', 'off');
                                    eval(['model.geom(''geom1'').feature(tmptag_mit1).selection(''input'').set({''',tmptag_mit0,''' ''ext1(6)''});']);
                                    model.geom('geom1').run(tmptag_mit1);
                                    tmp_getNDomains=model.geom('geom1').obj(tmptag_mit1).getNDomains;
                                    if tmp_getNDomains==0
                                        %                                     model.geom('geom1').feature(tmptag_mit0).active(false);
                                        %                                     model.geom('geom1').feature(tmptag_mit1).active(false);
                                        model.geom('geom1').feature.remove(tmptag_mit0);
                                        model.geom('geom1').feature.remove(tmptag_mit1);
                                    else
                                        model.geom('geom1').measure.selection.init(3);
                                        model.geom('geom1').measure.selection.set(tmptag_mit1, [1]);
                                        tmp_vol_mit=model.geom('geom1').measure().getVolume();
                                        if (vol_mit0/2-tmp_vol_mit)/(vol_mit0/2)<1e-6 % take this as equal
                                            %enable & create selection & record
                                            count_mit=count_mit+1;
                                            tmptag_mitall(count_mit)=cellstr(tmptag_mit1);
                                            tmptag_mitall_del(count_mit)=cellstr(tmptag_mit0);
                                        else
                                            %disable
                                            %                                         model.geom('geom1').feature(tmptag_mit0).active(false);
                                            %                                         model.geom('geom1').feature(tmptag_mit1).active(false);
                                            model.geom('geom1').feature.remove(tmptag_mit0);
                                            model.geom('geom1').feature.remove(tmptag_mit1);
                                        end
                                    end
                                end
                                %##operation int mito error
                                size(tmptag_mitall);
                                if size(tmptag_mitall,1)==0
                                    %                                 model.geom('geom1').feature(tmptag_vac0).active(false);
                                    %                                 model.geom('geom1').feature(tmptag_vac1).active(false);
                                    model.geom('geom1').feature.remove(tmptag_vac0);
                                    model.geom('geom1').feature.remove(tmptag_vac1);
                                    %                                 model.geom('geom1').feature(tmptag_chl0).active(false);
                                    %                                 model.geom('geom1').feature(tmptag_chl1).active(false);
                                    model.geom('geom1').feature.remove(tmptag_chl0);
                                    model.geom('geom1').feature.remove(tmptag_chl1);
                                else
                                    %%copy mov4(whole MSC boundary)
                                    tmptag_msw0=['copy',num2str(count_copy)];
                                    model.geom('geom1').feature.create(tmptag_msw0, 'Copy');
                                    count_copy=count_copy+1;
                                    eval(['model.geom(''geom1'').feature(tmptag_msw0).set(''displx'',''col_dx*',num2str(loop_j),'-floor(',num2str(loop_j),'/3)*row_dx+(',num2str(loop_i),'-1)*row_dx+scale_x*1*(',num2str(tmp_dx),')'');']);
                                    eval(['model.geom(''geom1'').feature(tmptag_msw0).set(''disply'',''col_dy*',num2str(loop_j),'-floor(',num2str(loop_j),'/3)*row_dy+(',num2str(loop_i),'-1)*row_dy'');']);
                                    model.geom('geom1').feature(tmptag_msw0).selection('input').set(tag_types{ms_type,4});
                                    %difference whole cell object-new vacuole&mito&chlo = new cyto [create selection]
                                    tmptag_msw1=['int',num2str(count_int)];
                                    model.geom('geom1').feature.create(tmptag_msw1, 'Intersection');
                                    count_int=count_int+1;
                                    model.geom('geom1').feature(tmptag_msw1).set('intbnd', 'off');
                                    model.geom('geom1').feature(tmptag_msw1).set('keep', 'on');
                                    clear tmptag_set
                                    tmptag_set(1)=cellstr(tmptag_msw0);
                                    tmptag_set(2)=cellstr('ext1(5)');
                                    model.geom('geom1').feature(tmptag_msw1).selection('input').set(tmptag_set);
                                    model.geom('geom1').run(tmptag_msw1);
                                    if model.geom('geom1').obj(tmptag_msw1).getNDomains==0
                                        %                                     model.geom('geom1').feature(tmptag_vac0).active(false);
                                        %                                     model.geom('geom1').feature(tmptag_vac1).active(false);
                                        model.geom('geom1').feature.remove(tmptag_vac0);
                                        model.geom('geom1').feature.remove(tmptag_vac1);
                                        %                                     model.geom('geom1').feature(tmptag_chl0).active(false);
                                        %                                     model.geom('geom1').feature(tmptag_chl1).active(false);
                                        model.geom('geom1').feature.remove(tmptag_chl0);
                                        model.geom('geom1').feature.remove(tmptag_chl1);
                                        for loop_k=1:size(tmptag_mitall,2)
                                            %                                         model.geom('geom1').feature(tmptag_mitall(loop_k)).active(false);
                                            %                                         model.geom('geom1').feature(tmptag_mitall_del(loop_k)).active(false);
                                            model.geom('geom1').feature.remove(tmptag_mitall(loop_k));
                                            model.geom('geom1').feature.remove(tmptag_mitall_del(loop_k));
                                        end
                                        %                                     model.geom('geom1').feature(tmptag_msw0).active(false);
                                        %                                     model.geom('geom1').feature(tmptag_msw1).active(false);
                                        model.geom('geom1').feature.remove(tmptag_msw0);
                                        model.geom('geom1').feature.remove(tmptag_msw1);
                                    else
                                        model.geom('geom1').measure.selection.init(3);
                                        model.geom('geom1').measure.selection.set(tmptag_msw1, [1]);
                                        tmp_vol1=model.geom('geom1').measure().getVolume();
                                        model.geom('geom1').measure.selection.init(3);
                                        model.geom('geom1').measure.selection.set(tmptag_msw0, [1]);
                                        tmp_vol0=model.geom('geom1').measure().getVolume();
                                        if tmp_vol1>tmp_vol0
                                            %                                         model.geom('geom1').feature(tmptag_vac0).active(false);
                                            %                                         model.geom('geom1').feature(tmptag_vac1).active(false);
                                            model.geom('geom1').feature.remove(tmptag_vac0);
                                            model.geom('geom1').feature.remove(tmptag_vac1);
                                            %                                         model.geom('geom1').feature(tmptag_chl0).active(false);
                                            %                                         model.geom('geom1').feature(tmptag_chl1).active(false);
                                            model.geom('geom1').feature.remove(tmptag_chl0);
                                            model.geom('geom1').feature.remove(tmptag_chl1);
                                            for loop_k=1:size(tmptag_mitall,2)
                                                %                                             model.geom('geom1').feature(tmptag_mitall(loop_k)).active(false);
                                                %                                             model.geom('geom1').feature(tmptag_mitall_del(loop_k)).active(false);
                                                model.geom('geom1').feature.remove(tmptag_mitall(loop_k));
                                                model.geom('geom1').feature.remove(tmptag_mitall_del(loop_k));
                                            end
                                            %                                         model.geom('geom1').feature(tmptag_msw0).active(false);
                                            %                                         model.geom('geom1').feature(tmptag_msw1).active(false);
                                            model.geom('geom1').feature.remove(tmptag_msw0);
                                            model.geom('geom1').feature.remove(tmptag_msw1);
                                        else
                                            tmptag_mscyt=['dif',num2str(count_dif)];
                                            model.geom('geom1').feature.create(tmptag_mscyt, 'Difference');
                                            count_dif=count_dif+1;
                                            model.geom('geom1').feature(tmptag_mscyt).set('keep', 'on');
                                            model.geom('geom1').feature(tmptag_mscyt).set('intbnd', 'off');
                                            model.geom('geom1').feature(tmptag_mscyt).selection('input').set(tmptag_msw1);
                                            clear tmptag_set
                                            tmptag_set=[tmptag_mitall,tmptag_vac1,tmptag_chl1];
                                            model.geom('geom1').feature(tmptag_mscyt).selection('input2').set(tmptag_set);
                                            %model.geom('geom1').feature(tmptag_mscyt).set('createselection', true);%create selection before run
                                            %##operation dif error
                                            model.geom('geom1').run(tmptag_mscyt);
                                            if model.geom('geom1').obj(tmptag_mscyt).getNDomains==0
                                                %                                             model.geom('geom1').feature(tmptag_vac0).active(false);
                                                %                                             model.geom('geom1').feature(tmptag_vac1).active(false);
                                                model.geom('geom1').feature.remove(tmptag_vac0);
                                                model.geom('geom1').feature.remove(tmptag_vac1);
                                                %                                             model.geom('geom1').feature(tmptag_chl0).active(false);
                                                %                                             model.geom('geom1').feature(tmptag_chl1).active(false);
                                                model.geom('geom1').feature.remove(tmptag_chl0);
                                                model.geom('geom1').feature.remove(tmptag_chl1);
                                                for loop_k=1:size(tmptag_mitall,2)
                                                    %                                                 model.geom('geom1').feature(tmptag_mitall(loop_k)).active(false);
                                                    %                                                 model.geom('geom1').feature(tmptag_mitall_del(loop_k)).active(false);
                                                    model.geom('geom1').feature.remove(tmptag_mitall(loop_k));
                                                    model.geom('geom1').feature.remove(tmptag_mitall_del(loop_k));
                                                end
                                                %                                             model.geom('geom1').feature(tmptag_msw0).active(false);
                                                %                                             model.geom('geom1').feature(tmptag_msw1).active(false);
                                                %                                             model.geom('geom1').feature(tmptag_mscyt).active(false);
                                                model.geom('geom1').feature.remove(tmptag_msw0);
                                                model.geom('geom1').feature.remove(tmptag_msw1);
                                                model.geom('geom1').feature.remove(tmptag_mscyt);
                                            else
                                                model.geom('geom1').measure.selection.init(3);
                                                model.geom('geom1').measure.selection.set(tmptag_mscyt, [1]);
                                                tmp_vol1=model.geom('geom1').measure().getVolume();
                                                model.geom('geom1').measure.selection.init(3);
                                                model.geom('geom1').measure.selection.set(tmptag_msw1, [1]);
                                                tmp_vol0=model.geom('geom1').measure().getVolume();
                                                if tmp_vol1>tmp_vol0
                                                    %                                                 model.geom('geom1').feature(tmptag_vac0).active(false);
                                                    %                                                 model.geom('geom1').feature(tmptag_vac1).active(false);
                                                    model.geom('geom1').feature.remove(tmptag_vac0);
                                                    model.geom('geom1').feature.remove(tmptag_vac1);
                                                    %                                                 model.geom('geom1').feature(tmptag_chl0).active(false);
                                                    %                                                 model.geom('geom1').feature(tmptag_chl1).active(false);
                                                    model.geom('geom1').feature.remove(tmptag_chl0);
                                                    model.geom('geom1').feature.remove(tmptag_chl1);
                                                    for loop_k=1:size(tmptag_mitall,2)
                                                        %                                                     model.geom('geom1').feature(tmptag_mitall(loop_k)).active(false);
                                                        %                                                     model.geom('geom1').feature(tmptag_mitall_del(loop_k)).active(false);
                                                        model.geom('geom1').feature.remove(tmptag_mitall(loop_k));
                                                        model.geom('geom1').feature.remove(tmptag_mitall_del(loop_k));
                                                    end
                                                    %                                                 model.geom('geom1').feature(tmptag_msw0).active(false);
                                                    %                                                 model.geom('geom1').feature(tmptag_msw1).active(false);
                                                    %                                                 model.geom('geom1').feature(tmptag_mscyt).active(false);
                                                    model.geom('geom1').feature.remove(tmptag_msw0);
                                                    model.geom('geom1').feature.remove(tmptag_msw1);
                                                    model.geom('geom1').feature.remove(tmptag_mscyt);
                                                else
                                                    %%union all mito
                                                    tmptag_uni=['uni',num2str(count_uni)];
                                                    model.geom('geom1').feature.create(tmptag_uni, 'Union');
                                                    count_uni=count_uni+1;
                                                    model.geom('geom1').feature(tmptag_uni).selection('input').set(tmptag_mitall);
                                                    %model.geom('geom1').feature(tmptag_uni).set('createselection', true);%save tag
                                                    count_list_cresel=count_list_cresel+1;
                                                    list_cresel(count_list_cresel)=cellstr(tmptag_uni);
                                                    tmptag_mituni=tmptag_uni;
                                                    prep_CreSel(count_list_cresel,1)=count_list_cresel;
                                                    prep_CreSel(count_list_cresel,2)=3;
                                                    %%%create selection and tmptags for intercellular air space
                                                    %model.geom('geom1').feature(tmptag_vac1).set('createselection', true);%save tag
                                                    count_list_cresel=count_list_cresel+1;
                                                    list_cresel(count_list_cresel)=cellstr(tmptag_vac1);
                                                    prep_CreSel(count_list_cresel,1)=count_list_cresel;
                                                    prep_CreSel(count_list_cresel,2)=4;
                                                    %model.geom('geom1').feature(tmptag_chl1).set('createselection', true);%save tag
                                                    count_list_cresel=count_list_cresel+1;
                                                    list_cresel(count_list_cresel)=cellstr(tmptag_chl1);
                                                    prep_CreSel(count_list_cresel,1)=count_list_cresel;
                                                    prep_CreSel(count_list_cresel,2)=2;
                                                    %model.geom('geom1').feature(tmptag_mscyt).set('createselection', true);%save tag
                                                    count_list_cresel=count_list_cresel+1;
                                                    list_cresel(count_list_cresel)=cellstr(tmptag_mscyt);
                                                    prep_CreSel(count_list_cresel,1)=count_list_cresel;
                                                    prep_CreSel(count_list_cresel,2)=1;
                                                    %%%change selection name
                                                    tmp_str=['ms',num2str(count_MSC),'_vac'];
                                                    model.geom('geom1').feature(tmptag_vac1).name(tmp_str);
                                                    tmp_str=['ms',num2str(count_MSC),'_chl'];
                                                    model.geom('geom1').feature(tmptag_chl1).name(tmp_str);
                                                    tmp_str=['ms',num2str(count_MSC),'_mits'];
                                                    model.geom('geom1').feature(tmptag_mituni).name(tmp_str);
                                                    tmp_str=['ms',num2str(count_MSC),'_cyt'];
                                                    model.geom('geom1').feature(tmptag_mscyt).name(tmp_str);
                                                    tmptag_IASinput2(count_MSC)=cellstr(tmptag_msw0);
                                                    count_MSC;
                                                    count_MSC=count_MSC+1;
                                                    %%%delete
                                                    tmptag_del=['del',num2str(count_del)];
                                                    model.geom('geom1').feature.create(tmptag_del, 'Delete');
                                                    count_del=count_del+1;
                                                    model.geom('geom1').feature(tmptag_del).selection('input').init;
                                                    clear tmptag_set
                                                    tmptag_set=[tmptag_vac0,tmptag_chl0,tmptag_mitall_del,tmptag_msw1];
                                                    model.geom('geom1').feature(tmptag_del).selection('input').set(tmptag_set);
                                                    %extend prep_CreSel
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%% save geom
%mphsave(model,'tmp_geom43.mph')

%%create IAS
tmptag_uni=['uni',num2str(count_uni)];
model.geom('geom1').feature.create(tmptag_uni, 'Union');
count_uni=count_uni+1;
model.geom('geom1').feature(tmptag_uni).set('intbnd', 'off');
model.geom('geom1').feature(tmptag_uni).selection('input').set(tmptag_IASinput2);

tmptag_IAS=['dif',num2str(count_dif)];
model.geom('geom1').feature.create(tmptag_IAS, 'Difference');
count_dif=count_dif+1;
model.geom('geom1').feature(tmptag_IAS).set('keep', 'on');
model.geom('geom1').feature(tmptag_IAS).set('intbnd', 'off');
model.geom('geom1').feature(tmptag_IAS).selection('input').set('ext1(5)');
model.geom('geom1').feature(tmptag_IAS).selection('input2').set(tmptag_uni);
%model.geom('geom1').feature(tmptag_IAS).set('repairtol', '1e-6');
model.geom('geom1').feature(tmptag_IAS).set('repairtoltype', 'auto');
% model.geom('geom1').feature(tmptag_IAS).set('createselection', true);%create selection before run
% tmp_str=['ias'];
% model.geom('geom1').feature(tmptag_IAS).name(tmp_str);
%##operation dif error
model.geom('geom1').run(tmptag_IAS);
if model.geom('geom1').obj(tmptag_IAS).getNDomains==0
    display('Error: IAS 1');
    %break;
else
    model.geom('geom1').measure.selection.init(3);
    model.geom('geom1').measure.selection.set(tmptag_IAS, [1]);
    tmp_vol1=model.geom('geom1').measure().getVolume();
    model.geom('geom1').measure.selection.init(3);
    model.geom('geom1').measure.selection.set('ext1(5)', [1]);
    tmp_vol0=model.geom('geom1').measure().getVolume();
    if tmp_vol1>tmp_vol0
        display('Error: IAS 2');
        %break;
    else
        tmptag_del=['del',num2str(count_del)];
        model.geom('geom1').feature.create(tmptag_del, 'Delete');
        count_del=count_del+1;
        model.geom('geom1').feature(tmptag_del).selection('input').init;
        clear tmptag_set
        %tmptag_set=['ext1(5)',tmptag_IASinput2];
        tmptag_set=['ext1(5)',cellstr(tmptag_uni)];
        model.geom('geom1').feature(tmptag_del).selection('input').set(tmptag_set);

        
        %%add stomata
        model.param.set('SMbaseU_width', '5e-6[m]');
        model.param.set('SMbaseU_depth', '1e-6[m]');
        model.geom('geom1').feature.create('blk1', 'Block');
        model.geom('geom1').feature('blk1').set('size', {'SMbaseU_width' 'SMbaseU_depth' 'LEAF_z'});
        model.geom('geom1').feature('blk1').set('pos', {'BSE_width/2+EPU_width/2-SMbaseU_width/2' 'EPL_thick+MST_thickatvein' '0'});
        model.param.set('SMU_width','3e-6[m]');
        model.param.set('SMU_depth','1e-6[m]');
        model.param.set('SMU_height','3e-6[m]');
        model.geom('geom1').feature.create('blk2', 'Block');
        model.geom('geom1').feature('blk2').set('size', {'SMU_width' 'SMU_depth' 'SMU_height'});
        model.geom('geom1').feature('blk2').set('pos', {'BSE_width/2+EPU_width/2-SMU_width/2' 'EPL_thick+MST_thickatvein+SMbaseU_depth' 'LEAF_z/2-SMU_height/2'});
        %model.geom('geom1').feature('blk2').set('createselection', true);%create selection before run
        count_list_cresel=count_list_cresel+1;
        list_cresel(count_list_cresel)=cellstr('blk2');
        prep_CreSel(count_list_cresel,1)=count_list_cresel;
        prep_CreSel(count_list_cresel,2)=6;
        tmp_str='smu';%upper stomata
        model.geom('geom1').feature('blk2').name(tmp_str);
        model.param.set('SMbaseL_width', '5e-6[m]');
        model.param.set('SMbaseL_depth', '1e-6[m]');
        model.geom('geom1').feature.create('blk3', 'Block');
        model.geom('geom1').feature('blk3').set('size', {'SMbaseL_width' 'SMbaseL_depth' 'LEAF_z'});
        model.geom('geom1').feature('blk3').set('pos', {'BSE_width/2+EPU_width/2+BUT_width/4-SMbaseL_width/2' '-SMbaseL_depth+EPL_thick' '0'});
        model.param.set('SML_width','3e-6[m]');
        model.param.set('SML_depth','1e-6[m]');
        model.param.set('SML_height','3e-6[m]');
        model.geom('geom1').feature.create('blk4', 'Block');
        model.geom('geom1').feature('blk4').set('size', {'SML_width' 'SML_depth' 'SML_height'});
        model.geom('geom1').feature('blk4').set('pos', {'BSE_width/2+EPU_width/2+BUT_width/4-SML_width/2' '-SMbaseL_depth-SML_depth+EPL_thick' 'LEAF_z/2-SML_height/2'});
        model.geom('geom1').run('blk4');
        %model.geom('geom1').feature('blk4').set('createselection', true);%create selection before run
        count_list_cresel=count_list_cresel+1;
        list_cresel(count_list_cresel)=cellstr('blk4');
        prep_CreSel(count_list_cresel,1)=count_list_cresel;
        prep_CreSel(count_list_cresel,2)=7;
        tmp_str='sml';%upper stomata
        model.geom('geom1').feature('blk4').name(tmp_str);
        tmptag_uni=['uni',num2str(count_uni)];
        model.geom('geom1').feature.create(tmptag_uni, 'Union');
        count_uni=count_uni+1;
        clear tmptag_set
        tmptag_set={tmptag_IAS,'blk1','blk3'};
        model.geom('geom1').feature(tmptag_uni).selection('input').set(tmptag_set);
        model.geom('geom1').feature(tmptag_uni).set('intbnd', 'off');
        %model.geom('geom1').feature(tmptag_uni).set('createselection', true);%create selection before run
        tmp_str='ias';
        model.geom('geom1').feature(tmptag_uni).name(tmp_str);
        
        count_list_cresel=count_list_cresel+1;
        list_cresel(count_list_cresel)=cellstr(tmptag_uni);
        prep_CreSel(count_list_cresel,1)=count_list_cresel;
        prep_CreSel(count_list_cresel,2)=5;
    end
end
%%delete mov1 mov2 mov3 mov4
tmptag_del=['del',num2str(count_del)];
model.geom('geom1').feature.create(tmptag_del, 'Delete');
count_del=count_del+1;
model.geom('geom1').feature(tmptag_del).selection('input').init;
tmp_test=reshape(tag_types,1,[]);
model.geom('geom1').feature(tmptag_del).selection('input').set(tmp_test);
model.geom('geom1').run(tmptag_del);

tmptag_del=['del',num2str(count_del)];
model.geom('geom1').feature.create(tmptag_del, 'Delete');
count_del=count_del+1;
model.geom('geom1').feature(tmptag_del).selection('input').init;
model.geom('geom1').feature(tmptag_del).selection('input').set({'ext1(1)' 'ext1(3)' 'ext1(4)' 'ext1(6)'});
model.geom('geom1').run(tmptag_del);


mphsave(model,['tmpCK_geom_cad_nocresel_',num2str(loop_idx),'_',num2str(loop_idy),'.mph'])
%toc;
%ModelUtil.disconnect;
%exit;

