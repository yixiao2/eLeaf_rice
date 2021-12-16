% eLeaf: 3D model of rice leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 1.2.6

clear all;%tic
import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model_IAS2D');
load parainput.mat
%%%%%%%%%%debug%%%%%%%%%%%%%%
%IAS_2D_input=1;
% cell_length=14.7*3;
% cell_height=10.6*3;
%%%%%%%%%%%%%%%

model.param.set('BSE_width', [num2str(BSE_width,GLB_digits),'[m]']);
model.param.set('EPL_thick', [num2str(EPL_thick,GLB_digits),'[m]'], 'thickness of lower epidermis');
model.param.set('MST_thickatvein', [num2str(MST_thickatvein,GLB_digits),'[m]'], 'mesophyll thickness at vein');
model.param.set('EPU_width', [num2str(EPU_width,GLB_digits),'[m]']);
model.param.set('BUT_width', [num2str(BUT_width,GLB_digits),'[m]']);
model.param.set('MST_thickatBU', [num2str(MST_thickatBU,GLB_digits),'[m]']);
model.param.set('VEIN_r', [num2str(VEIN_r,GLB_digits),'[m]']);
model.param.set('BS_thick', [num2str(BS_thick,GLB_digits),'[m]']);
model.param.set('VEIN_l', [num2str(VEIN_l,GLB_digits),'[m]']);
model.param.set('EPU_thick', [num2str(EPU_thick,GLB_digits),'[m]']);
model.param.set('BU_thick', [num2str(BU_thick,GLB_digits),'[m]']);
model.param.set('MS_dthick', [num2str(MS_dthick,GLB_digits),'[m]']);
model.param.set('BS_dthick', [num2str(BS_dthick,GLB_digits),'[m]']);
model.param.set('BS_mito_l', [num2str(BS_mito_l,GLB_digits),'[m]']);
model.param.set('BS_mito_r', [num2str(BS_mito_r,GLB_digits),'[m]']);
model.modelNode.create('mod1');

model.geom.create('geom1', 3);
%model.geom('geom1').geomRep('comsol');
%% DEFAULT relative repair tolerance
%model.geom('geom1').repairTol(1.0E-8);
model.geom('geom1').repairTolType('auto');
model.geom('geom1').feature.create('wp1', 'WorkPlane');
model.geom('geom1').feature('wp1').geom.feature.create('b1', 'BezierPolygon');
model.geom('geom1').feature('wp1').geom.feature.create('pol1', 'Polygon');
model.geom('geom1').feature('wp1').geom.feature.create('c1', 'Circle');
model.geom('geom1').feature('wp1').geom.feature.create('b2', 'BezierPolygon');
model.geom('geom1').feature('wp1').geom.feature.create('dif2', 'Difference');
model.geom('geom1').feature('wp1').geom.feature.create('dif4', 'Difference');
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
model.geom('geom1').feature('wp1').geom.feature('dif2').name('BSE');
model.geom('geom1').feature('wp1').geom.feature('dif2').set('keep', true);
model.geom('geom1').feature('wp1').geom.feature('dif2').set('intbnd', false);
model.geom('geom1').feature('wp1').geom.feature('dif2').set('edge', 'all');
model.geom('geom1').feature('wp1').geom.feature('dif2').selection('input2').set({'c1'});
model.geom('geom1').feature('wp1').geom.feature('dif2').selection('input').set({'pol1'});
model.geom('geom1').feature('wp1').geom.feature('dif4').name('MS');
model.geom('geom1').feature('wp1').geom.feature('dif4').set('keep', true);
model.geom('geom1').feature('wp1').geom.feature('dif4').set('intbnd', false);
model.geom('geom1').feature('wp1').geom.feature('dif4').set('edge', 'all');
model.geom('geom1').feature('wp1').geom.feature('dif4').selection('input2').set({'c1'});
model.geom('geom1').feature('wp1').geom.feature('dif4').selection('input').set({'b1'});
model.geom('geom1').feature('wp1').geom.feature('del1').selection('input').init;
% model.geom('geom1').feature('wp1').geom.feature('del1').selection('input').set({'b1' 'b2' 'c1' 'pol1'});
model.geom('geom1').feature('wp1').geom.feature('del1').selection('input').set({'b1' 'pol1'});

%% fill in outline of mesophyll cell
%Nlobe=8;% number of lobes; optional 4+2n (n>=0)
model.param.set('Nlobe',num2str(Nlobe));
Tlobe=Nlobe/2;% type of lobe
model.param.set('Tlobe',num2str(Tlobe));
%dw=0.05;
model.param.set('dw', num2str(dw));
dse=0.5;
model.param.set('dse', num2str(dse));
%dvaz=0.1;
model.param.set('dvaz', num2str(dvaz));
%cell_length=14.7;%unit um
%cell_height=10.6;%unit um
%cell_volume=1783;%unit um^3
%calculate cell_thick in z axis
cell_thick=round(cell_volume*3/4/pi/cell_length*2/cell_height*2*2,,GLB_digits,'significant');
model.param.set('cell_length',[num2str(cell_length,GLB_digits),'[um]']);
model.param.set('cell_height',[num2str(cell_height,GLB_digits),'[um]']);
model.param.set('cell_thick',[num2str(cell_thick,GLB_digits),'[um]']);
%mit_r=0.07;
%mit_l=0.05;
model.param.set('mit_r',num2str(mit_r));
model.param.set('mit_l',num2str(mit_l));
model.param.set('lmbd','0.2');  
model.param.set('rho','0');
%epsln=0.2;
if Nlobe==4
    model.param.set('epsln','0');
else
    model.param.set('epsln',num2str(epsln));
end
model.param.set('scale_x','(cell_length/2)/(1/2*Tlobe)');
model.param.set('scale_y','(cell_height/2)/(1/(2*sqrt(3))*(Tlobe+1))');
tmp_Npts=14+4*(Tlobe-2)+2*Nlobe;
% dw=0.05;dse=0.5;
[tmp_wallString,tmp_chloString,tmp_chliString,tmp_IDXxse,tmp_IDXyse]=msoutline_v1_2(Nlobe,dw,dse,epsln);


% model.geom('geom1').feature.create('wp2', 'WorkPlane');
count_wppol=2;tmptag_wppol=['pol',num2str(count_wppol)];
model.geom('geom1').feature('wp1').geom.feature.create(tmptag_wppol, 'Polygon');
count_wppol=count_wppol+1;
model.geom('geom1').feature('wp1').geom.feature(tmptag_wppol).set('source', 'table');
eval(['model.geom(''geom1'').feature(''wp1'').geom.feature(tmptag_wppol).set(''table'', {',tmp_wallString,'});']);

%% distribute ms cell
MSC_length=cell_length*1e-6;MSC_height=cell_height*1e-6;
N_ms_cols=floor((EPU_width+BUT_width/2)/MSC_length)+2;
N_ms_rows=floor(MST_thickatvein/MSC_height)+2;

% random the initial point
Origin_x=BSE_width/2-MSC_length/2*rand;
Origin_y=EPL_thick-MSC_height/2*rand;
Origin_x=str2num(num2str(Origin_x,GLB_digits));
Origin_y=str2num(num2str(Origin_y,GLB_digits));

model.geom('geom1').feature('wp1').geom.create('mov1', 'Move');
model.geom('geom1').feature('wp1').geom.feature('mov1').selection('input').set(tmptag_wppol);
model.geom('geom1').feature('wp1').geom.feature('mov1').set('displx', num2str(Origin_x,GLB_digits));
model.geom('geom1').feature('wp1').geom.feature('mov1').set('disply', num2str(Origin_y,GLB_digits));
model.geom('geom1').feature('wp1').geom.run('mov1');
model.geom('geom1').feature('wp1').geom.run('mov1');

model.param.set('col_dx', '(Nlobe/2-1/2)*scale_x');
model.param.set('col_dy', '3/2/sqrt(3)*scale_y');
model.param.set('row_dx', '1/2*scale_x');
model.param.set('row_dy', '9/2/sqrt(3)*scale_y');

model.geom('geom1').feature('wp1').geom.measure.selection.init(2);
model.geom('geom1').feature('wp1').geom.measure.selection.set('mov1', [1]);
sa_ms=model.geom('geom1').feature('wp1').geom.measure().getVolume();
model.geom('geom1').feature('wp1').geom.measure.selection.init(2);
model.geom('geom1').feature('wp1').geom.measure.selection.set('b2', [1]);
sa_leaf=model.geom('geom1').feature('wp1').geom.measure().getVolume();
model.geom('geom1').feature('wp1').geom.measure.selection.init(2);
model.geom('geom1').feature('wp1').geom.measure.selection.set('dif4', [1]);
sa_msregion=model.geom('geom1').feature('wp1').geom.measure().getVolume();

sa_thres=0.45;%threshold for ms too small
count_copy=1;
count_int=1;
count_ms=1;
ms_distribution=zeros(N_ms_rows,N_ms_cols+1,4);
sum_sa_ms=0;
for loop_i=1:N_ms_rows
    for loop_j=0:N_ms_cols
        tmptag_copy=['copy',num2str(count_copy)];
        model.geom('geom1').feature('wp1').geom.create(tmptag_copy, 'Copy');
        model.geom('geom1').feature('wp1').geom.feature(tmptag_copy).selection('input').set({'mov1'});
        eval(['model.geom(''geom1'').feature(''wp1'').geom.feature(tmptag_copy).set(''displx'',''col_dx*',num2str(loop_j),'-floor(',num2str(loop_j),'/3)*row_dx+(',num2str(loop_i),'-1)*row_dx'');']);
        eval(['model.geom(''geom1'').feature(''wp1'').geom.feature(tmptag_copy).set(''disply'',''col_dy*',num2str(loop_j),'-floor(',num2str(loop_j),'/3)*row_dy+(',num2str(loop_i),'-1)*row_dy'');']);
        count_copy=count_copy+1;        
        model.geom('geom1').feature('wp1').geom.run(tmptag_copy);
        
        tmptag_int=['int',num2str(count_int)];        
        model.geom('geom1').feature('wp1').geom.feature.create(tmptag_int, 'Intersection');
        count_int=count_int+1;
        model.geom('geom1').feature('wp1').geom.feature(tmptag_int).set('keep', 'on');
        model.geom('geom1').feature('wp1').geom.feature(tmptag_int).set('intbnd', 'off');
        eval(['model.geom(''geom1'').feature(''wp1'').geom.feature(tmptag_int).selection(''input'').set({''',tmptag_copy,''' ''dif4''});']);
        model.geom('geom1').feature('wp1').geom.run(tmptag_int);
        
        tmp_getNDomains=model.geom('geom1').feature('wp1').geom.obj(tmptag_int).getNDomains;
        if tmp_getNDomains==0
            model.geom('geom1').feature('wp1').geom.feature.remove(tmptag_copy);
            model.geom('geom1').feature('wp1').geom.feature.remove(tmptag_int);
        else
            model.geom('geom1').feature('wp1').geom.measure.selection.init(2);
            model.geom('geom1').feature('wp1').geom.measure.selection.set(tmptag_int, [1]);
            tmpsa_int=model.geom('geom1').feature('wp1').geom.measure().getVolume();
            if tmpsa_int/sa_ms<sa_thres
                model.geom('geom1').feature('wp1').geom.feature.remove(tmptag_copy);
                model.geom('geom1').feature('wp1').geom.feature.remove(tmptag_int);
            else
                tmptag_msall(count_ms)=cellstr(tmptag_copy);
                count_ms=count_ms+1;
                if tmpsa_int/sa_ms>0.999
                    ms_distribution(loop_i,loop_j+1,1)=2;
                else
                    ms_distribution(loop_i,loop_j+1,1)=1;
                end
                sum_sa_ms=sum_sa_ms+tmpsa_int;
            end 
        end
    end
end
model.geom('geom1').feature('wp1').geom.create('del2', 'Delete');
model.geom('geom1').feature('wp1').geom.feature('del2').selection('input').init;
model.geom('geom1').feature('wp1').geom.feature('del2').selection('input').set(['mov1',tmptag_msall]);
%model.geom('geom1').feature('wp1').geom.run('del2');

IAS_2D_min_1=1-sum_sa_ms/sa_msregion;%IAS per ms_region
IAS_2D_max_1=IAS_2D_min_1+sa_ms/sa_msregion*sum(sum(ms_distribution(:,:,1)==2));
IAS_2D_min_2= (sa_msregion-sum_sa_ms)/sa_leaf;%IAS per leaf cross-section

    Nhole_IDX=find(ms_distribution(:,:,1)==2);
    tmp_page=-1*ones(size(ms_distribution(:,:,1)));
    tmp_page2=tmp_page;
    tmp_page(Nhole_IDX)=0;
    tmp_page2(Nhole_IDX)=0;
%     display(['Error: Exact IAD_2D_input cannot be reached since it is too small. The modeled IAS_2D is ',num2str(round(IAS_2D_min_1,GLB_digits,'significant')),' .']);    
    lmbd=0.25;rho=0;

ms_distribution(:,:,2)=tmp_page;
ms_distribution(:,:,3)=tmp_page2;

save tmp_IAS2D.mat -regexp '^(?!(model|ans)$).'
display('IAS2D matched.')
%toc
