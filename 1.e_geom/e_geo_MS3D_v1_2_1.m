% eLeaf: 3D model of rice leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 1.2.6

clear all;%clc

load parainput.mat

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model_3Dcell2');

model.name('automatic_cell_RD_COMSOL.mph');

%Nlobe=4;% number of lobes; optional 4+2n (n>=0)
model.param.set('Nlobe',num2str(Nlobe));
Tlobe=Nlobe/2;% type of lobe
model.param.set('Tlobe',num2str(Tlobe));
%dw=0.05;
model.param.set('dw', num2str(dw));
dse=0.5;
model.param.set('dse', num2str(dse));
model.param.set('dvaz', '0.05');
%cell_length=14.7;%unit um
%cell_height=10.6;%unit um
%cell_volume=1783;%unit um^3
%calculate cell_thick in z axis
cell_thick=round(cell_volume*3/4/pi/cell_length*2/cell_height*2*2,GLB_digits,'significant');
model.param.set('cell_length',[num2str(cell_length,GLB_digits),'[um]']);
model.param.set('cell_height',[num2str(cell_height,GLB_digits),'[um]']);
model.param.set('cell_thick',[num2str(cell_thick,GLB_digits),'[um]']);
%mit_r=0.07;
%mit_l=0.05;
model.param.set('mit_r',num2str(mit_r));
model.param.set('mit_l',num2str(mit_l));
model.param.set('lmbd',num2str(lmbd));
model.param.set('rho',num2str(rho));
%epsln=0.2;
if Nlobe==4
    model.param.set('epsln','0');
else
    model.param.set('epsln',num2str(epsln));
end

model.modelNode.create('mod1');

%% geometry
model.param.set('scale_x','(cell_length/2)/(1/2*Tlobe)');
scale_x=(cell_length/2)/(1/2*Tlobe);
model.param.set('scale_y','(cell_height/2)/(1/(2*sqrt(3))*(Tlobe+1))');
scale_y=(cell_height/2)/(1/(2*sqrt(3))*(Tlobe+1));
tmp_Npts=14+4*(Tlobe-2)+2*Nlobe;
[tmp_wallString,tmp_chloString,tmp_chliString,tmp_IDXxse,tmp_IDXyse]=msoutline_v1_2(Nlobe,dw,dse,epsln);

model.geom.create('geom1', 3);
model.geom('geom1').geomRep('comsol');
model.geom('geom1').feature.create('wp1', 'WorkPlane');
model.geom('geom1').feature.create('ext1', 'Extrude');
model.geom('geom1').feature.create('elp1', 'Ellipsoid');
model.geom('geom1').feature.create('elp2', 'Ellipsoid');
model.geom('geom1').feature.create('elp3', 'Ellipsoid');
model.geom('geom1').feature.create('elp4', 'Ellipsoid');
model.geom('geom1').feature.create('int1', 'Intersection');
model.geom('geom1').feature('wp1').geom.feature.create('pol1', 'Polygon');
model.geom('geom1').feature('wp1').geom.feature.create('pol2', 'Polygon');
model.geom('geom1').feature('wp1').geom.feature.create('pol3', 'Polygon');
model.geom('geom1').feature('wp1').set('quickz', '-cell_thick');
model.geom('geom1').feature('wp1').geom.feature('pol1').set('source', 'table');
eval(['model.geom(''geom1'').feature(''wp1'').geom.feature(''pol1'').set(''table'', {',tmp_wallString,'});']);
model.geom('geom1').feature('wp1').geom.feature('pol2').set('source', 'table');
eval(['model.geom(''geom1'').feature(''wp1'').geom.feature(''pol2'').set(''table'', {',tmp_chloString,'});']);
model.geom('geom1').feature('wp1').geom.feature('pol3').set('source', 'table');
eval(['model.geom(''geom1'').feature(''wp1'').geom.feature(''pol3'').set(''table'', {',tmp_chliString,'});']);
model.geom('geom1').feature('ext1').setIndex('distance', 'cell_thick*2', 0);
model.geom('geom1').feature('ext1').selection('input').set({'wp1'});
rab=1;
%flatness=1.5;
ea=rab*sqrt((Tlobe/2)^2/rab^2+25/4/3)*flatness;
eb=sqrt((Tlobe/2)^2/rab^2+25/4/3)*flatness;
ea=str2num(num2str(ea,GLB_digits));
eb=str2num(num2str(eb,GLB_digits));
eval(['model.geom(''geom1'').feature(''elp1'').set(''semiaxes'', {''scale_x*',num2str(ea,GLB_digits),''' ''scale_y*',num2str(eb,GLB_digits),''' ''cell_thick/2''});'])
eval(['model.geom(''geom1'').feature(''elp2'').set(''semiaxes'', {''scale_x*(',num2str(ea,GLB_digits),'-dw)'' ''scale_y*(',num2str(eb,GLB_digits),'-dw)'' ''cell_thick/2-dw*(scale_x+scale_y)/2''});'])
eval(['model.geom(''geom1'').feature(''elp3'').set(''semiaxes'', {''scale_x*(',num2str(ea,GLB_digits),'-dse*(1-dw+epsln/2))'' ''scale_y*(',num2str(eb,GLB_digits),'-dse*(1-dw+epsln/2))'' ''cell_thick/2-dse*(1-dw+epsln/2)*(scale_x+scale_y)/2''});'])
% eval(['model.geom(''geom1'').feature(''elp1'').set(''semiaxes'', {''1.7*scale_x*',num2str(Tlobe/3),''' ''1.7*scale_y*',num2str(max(1,Tlobe/3)),''' ''cell_thick/2''});'])
% eval(['model.geom(''geom1'').feature(''elp2'').set(''semiaxes'', {''1.7*scale_x*(',num2str(Tlobe/3),'-dw)'' ''1.7*scale_y*(',num2str(max(1,Tlobe/3)),'-dw)'' ''cell_thick/2*(1-dw)''});'])
if Nlobe==4
    eval(['model.geom(''geom1'').feature(''elp4'').set(''semiaxes'', {''-scale_x*(',tmp_IDXxse{1},'+(mit_r+mit_l)*2+mit_r+0.05)'' ''scale_y*(',tmp_IDXyse{3},'-(mit_r+mit_l)*2-0.05)'' ''cell_thick/2-(dse*(1-dw+epsln/2)+dvaz)*(scale_x+scale_y)/2''});'])%x: pt6; y:pt8
else
    eval(['model.geom(''geom1'').feature(''elp4'').set(''semiaxes'', {''-scale_x*(',tmp_IDXxse{6},'+mit_r+0.05)'' ''scale_y*(',tmp_IDXyse{8},'-(mit_r+mit_l)*2-0.05)'' ''cell_thick/2-(dse*(1-dw+epsln/2)+dvaz)*(scale_x+scale_y)/2''});'])%x: pt6; y:pt8
end
% if Nlobe==4
%     model.geom('geom1').feature('elp4').active(false);
% end
model.geom('geom1').feature('int1').set('face', 'all');
model.geom('geom1').feature('int1').set('intbnd', false);
model.geom('geom1').feature('int1').selection('input').set({'elp1' 'ext1(1)'});
%% match the measured cell volume
model.geom('geom1').run('int1');
model.geom('geom1').measure.selection.init(3);
model.geom('geom1').measure.selection.set('int1', [1]);
tmp_vol1=model.geom('geom1').measure().getVolume();
cell_thick
cell_thick=cell_thick*cell_volume/(tmp_vol1*1e18)
cell_thick=str2num(num2str(cell_thick,GLB_digits));
model.param.set('cell_thick',[num2str(cell_thick,GLB_digits),'[um]']);%end match cell volume

model.geom('geom1').feature.create('int2', 'Intersection');
model.geom('geom1').feature.create('int3', 'Intersection');
model.geom('geom1').feature('int2').set('face', 'all');
model.geom('geom1').feature('int2').set('intbnd', false);
model.geom('geom1').feature('int2').selection('input').set({'elp2' 'ext1(2)'});
model.geom('geom1').feature('int3').set('face', 'all');
model.geom('geom1').feature('int3').set('intbnd', false);
model.geom('geom1').feature('int3').selection('input').set({'elp3' 'ext1(3)'});

%% match the measured plastid volume
%plastid_volume=0.474;% unit is percentage
%tmp_vol1=cell_volume;
flag_ub=0;
flag_lb=0;
model.geom('geom1').run('int3');
model.geom('geom1').measure.selection.init(3);
model.geom('geom1').measure.selection.set('int1', [1]);
tmp_vol1=model.geom('geom1').measure().getVolume()*1e18
cell_volume
model.geom('geom1').measure.selection.init(3);
model.geom('geom1').measure.selection.set('int2', [1]);
tmp_vol2=model.geom('geom1').measure().getVolume()*1e18;
model.geom('geom1').measure.selection.init(3);
model.geom('geom1').measure.selection.set('int3', [1]);
tmp_vol3=model.geom('geom1').measure().getVolume()*1e18;
tmp_plastid_vol=(tmp_vol2-tmp_vol3)/tmp_vol1
tmp_new_dse=dse;

while abs((tmp_plastid_vol-plastid_volume)/plastid_volume)>0.01&&flag_ub<2&&flag_lb<2
    tmp_new_dse=tmp_new_dse+plastid_volume-tmp_plastid_vol;
    if tmp_new_dse>=0.8
        tmp_new_dse=0.8;
        flag_ub=flag_ub+1;
    elseif tmp_new_dse<=0.05
        tmp_new_dse=0.05;
        flag_lb=flag_lb+1;
    else
        flag_ub=0;
        flag_lb=0;
    end
    tmp_new_dse=str2num(num2str(tmp_new_dse,GLB_digits));
    model.param.set('dse',num2str(tmp_new_dse,GLB_digits));
    [tmp_wallString,tmp_chloString,tmp_chliString,tmp_IDXxse,tmp_IDXyse]=msoutline_v1_2(Nlobe,dw,tmp_new_dse,epsln);
    eval(['model.geom(''geom1'').feature(''wp1'').geom.feature(''pol1'').set(''table'', {',tmp_wallString,'});']);
    eval(['model.geom(''geom1'').feature(''wp1'').geom.feature(''pol2'').set(''table'', {',tmp_chloString,'});']);
    eval(['model.geom(''geom1'').feature(''wp1'').geom.feature(''pol3'').set(''table'', {',tmp_chliString,'});']);
    model.geom('geom1').runPre('int3');
    model.geom('geom1').feature('int3').selection('input').set({'elp3' 'ext1(3)'});
    model.geom('geom1').run('int3');
    model.geom('geom1').measure.selection.init(3);
    model.geom('geom1').measure.selection.set('int2', [1]);
    tmp_vol2=model.geom('geom1').measure().getVolume()*1e18;
    model.geom('geom1').measure.selection.init(3);
    model.geom('geom1').measure.selection.set('int3', [1]);
    tmp_vol3=model.geom('geom1').measure().getVolume()*1e18;
    tmp_plastid_vol=(tmp_vol2-tmp_vol3)/tmp_vol1;
    tmp_new_dse
    tmp_plastid_vol
end
dse=tmp_new_dse;
if flag_ub==2
    display(['The plastid volume set is too large. The modeled plastid volume is ',num2str(round(tmp_plastid_vol*100,GLB_digits,'significant')),'% of the cell volume.']);
end
if flag_lb==2
    display(['The plastid volume set is too small. The modeled plastid volume is ',num2str(round(tmp_plastid_vol*100,GLB_digits,'significant')),'% of the cell volume.']);
end
%end match plastid volume

model.param.set('col_dx', '(Nlobe/2-1/2)*scale_x');
model.param.set('col_dy', '3/2/sqrt(3)*scale_y');
model.param.set('row_dx', '1/2*scale_x');
model.param.set('row_dy', '9/2/sqrt(3)*scale_y');

%% Sm determines lmbd
%delete all objects except int1
model.geom('geom1').create('del1', 'Delete');
model.geom('geom1').feature('del1').selection('input').init;
model.geom('geom1').feature('del1').selection('input').set({'int2' 'int3' 'elp4'});
model.geom('geom1').run('del1');
%cylinder selection
model.geom('geom1').run('fin');
model.selection.create('cyl1', 'Cylinder');
model.selection('cyl1').set('entitydim', '2');

model.geom('geom1').measure.selection.init(2);
tmp_bnd_num=mphgetadj(model,'geom1','boundary','domain',1);
model.geom('geom1').measure.selection.set('fin', tmp_bnd_num);
S=model.geom('geom1').measure().getVolume()
tmp_sel=mphgetselection(model.selection('cyl1'));
tmp_bnd_num=tmp_sel.entities;
model.geom('geom1').measure.selection.init(2);
model.geom('geom1').measure.selection.set('fin', tmp_bnd_num);
min_CW_touching_air=model.geom('geom1').measure().getVolume()

%rho=0;
%SmS_max=1-(0.4014);%require manual input
SmS_max=1-(min_CW_touching_air/S);%require manual input
if Nlobe==4
    tmp_ratio1=18/4;
    tmp_ratio2=10/4;
else
    tmp_ratio1=(24+(Tlobe-3)*8)/(18+(Tlobe-3)*4);
    tmp_ratio2=12/(18+(Tlobe-3)*4);
end
lmbd=(1-SmS/SmS_max)/(tmp_ratio1);
if lmbd>0.5
    lmbd=1-SmS/SmS_max/(tmp_ratio2);
    display('lmbd>0.5');
end

%% porosity determines leaf depth
%IAS_3D_input=0.072
%model.geom('geom1').feature('wp1').geom.measure.selection.init(2);
%model.geom('geom1').feature('wp1').geom.measure.selection.set('pol1', [1]);
%tmp_area=model.geom('geom1').feature('wp1').geom.measure().getVolume();
tmp_area=sqrt(3)/2*(4+3*(Tlobe-2))*scale_x*scale_y*1e-12;
leaf_z0=tmp_vol1*1e-18/(1-IAS_3D_input)/tmp_area/2;
cell_thick/2*1e-6;
%adjust rho to make the values below comparable??
(leaf_z0-cell_thick/2*1e-6)*1e6;
1/sqrt(3)*rho*(scale_x+scale_y)/2;

save tmp_MS3D.mat -regexp '^(?!(model|ans)$).'
display('MS3D matched.')

%toc
