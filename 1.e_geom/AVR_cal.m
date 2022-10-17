% eLeaf: 3D model of rice leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 1.2.6

%% calculate AVR
%% prepare for ray tracing
%% modified from e_physics
model=mphload('tmpCK_geomIP_cad_mesh_cresel.mph');
load save_e_geom.mat
%% pre-process to identity pair; classify all the pairs to 1)wall-IAS;
%% 2)wall-wall; 3)cyto_o-chlo; 4)cyto_i-chlo; 5)cyto_i-mit; 6)cyto_i-vac.
%% *IAS is short for Intercellular Air Space
NDom=model.geom('geom1').getNDomains();% number of domains
NBnd=model.geom('geom1').getNBoundaries();% number of boundaries

if mod(model.selection.size(),4)==0
    NCreSel=model.selection.size()/4;% number of create selections
    TagsCreSel=model.selection.tags();
else
    Display('Error in create selection.');
    %break;
end
NPairs=model.pair.size();% number of identity pairs
TagsPairs=model.pair.tags();

prep_dom=ones(NDom,3)*(-1);% initialize the matrix to store the mapping
prep_bnd=ones(NBnd,4)*(-1);

%% Question: how to generate the following prep_CreSel automatically
%prep_CreSel=[1,4;2,1;3,2;4,1;5,3];%first column 1:5=1:NCreSel;second column 1=cyto;2=chlo;3=mito;4=vacu;
prep_pairs(:,1)=1:NPairs;

%% generate prep_dom matrix
prep_dom(:,1)=1:NDom;
count_chl=0;
count_mit=0;
for tmp_loop=1:NCreSel
    tmp_char=char(TagsCreSel(tmp_loop*4));% *4 is because each created selection corresponds to 4 selections (pts,edg,bnd,dom)
    tmp_domSet=eval(['model.selection(''',tmp_char,''').inputEntities']);
    prep_dom(tmp_domSet,2)=tmp_loop;%for tmp_loop2=1:size(tmp_domSet,1);Then prep_dom+prep_CreSel-->biochemical&physical parameters
    prep_dom(tmp_domSet,3)=prep_CreSel(tmp_loop,2);
    if prep_CreSel(tmp_loop,2)==2
        count_chl=count_chl+1;
        tmp_str=['mschlo_',num2str(count_chl)];
        model.selection.create(tmp_str, 'Explicit');
        model.selection(tmp_str).geom('geom1', 3);
        model.selection(tmp_str).set(tmp_domSet);
        model.selection(tmp_str).name(tmp_str);
    end
    if prep_CreSel(tmp_loop,2)==3
        count_mit=count_mit+1;
        tmp_str=['msmito_',num2str(count_mit)];
        model.selection.create(tmp_str, 'Explicit');
        model.selection(tmp_str).geom('geom1', 3);
        model.selection(tmp_str).set(tmp_domSet);
        model.selection(tmp_str).name(tmp_str);
    end
end
if(size(find(prep_dom==-1),1)~=0)
    disp('Error: certain DOMAIN is not included in the Create Selection')
    %break;
end

%% generate prep_bnd matrix
prep_bnd(:,1)=1:NBnd;

for tmp_loop=1:NCreSel
    tmp_char=char(TagsCreSel(tmp_loop*4-1));% *4-1 is because each created selection corresponds to 4 selections (pts,edg,bnd,dom)
    tmp_bndSet=eval(['model.selection(''',tmp_char,''').inputEntities']);
    prep_bnd(tmp_bndSet,2)=tmp_loop;
end
if(size(find(prep_bnd(:,2)==-1),1)~=0)
    display('Error: certain BOUNDARY is not included in the Create Selection')
    %break;
end

CreSelType2pairtype=zeros(10,10);
CreSelType2pairtype(1,1)=1;
CreSelType2pairtype(1,2)=2;CreSelType2pairtype(2,1)=2;
CreSelType2pairtype(1,3)=3;CreSelType2pairtype(3,1)=3;
CreSelType2pairtype(1,4)=4;CreSelType2pairtype(4,1)=4;
CreSelType2pairtype(1,5)=5;CreSelType2pairtype(5,1)=5;
CreSelType2pairtype(5,6)=6;CreSelType2pairtype(6,5)=6;
CreSelType2pairtype(5,7)=7;CreSelType2pairtype(7,5)=7;
CreSelType2pairtype(5,8)=8;CreSelType2pairtype(8,5)=8;
CreSelType2pairtype(1,8)=9;CreSelType2pairtype(8,1)=9;
CreSelType2pairtype(8,9)=10;CreSelType2pairtype(9,8)=10;
CreSelType2pairtype(8,10)=11;CreSelType2pairtype(10,8)=11;

for tmp_loop=1:NPairs
    tmp_char=char(TagsPairs(tmp_loop));
    tmp_pairSrc=eval(['model.pair(''',tmp_char,''').source().entities(2)']);
    tmp_pairDst=eval(['model.pair(''',tmp_char,''').destination().entities(2)']);
    prep_bnd(tmp_pairSrc,3)=tmp_loop;
    prep_bnd(tmp_pairDst,3)=tmp_loop;
    tmp_IDXdom1=mphgetadj(model,'geom1','domain','boundary',tmp_pairSrc(1));
    tmp_IDXdom2=mphgetadj(model,'geom1','domain','boundary',tmp_pairDst(1));
    tmp_IDXCreSelType1=prep_dom(tmp_IDXdom1,3);
    tmp_IDXCreSelType2=prep_dom(tmp_IDXdom2,3);

    tmp_pairtype=CreSelType2pairtype(tmp_IDXCreSelType1,tmp_IDXCreSelType2);
    prep_bnd(tmp_pairSrc,4)=tmp_pairtype;
    prep_bnd(tmp_pairDst,4)=tmp_pairtype;
    prep_pairs(tmp_loop,2)=tmp_pairtype;

    % prep_bnd(:,5) is subtype of column 4, specifically for
    % wall_ms-wall_ms VS wall_bs-wall_ms. This will be activated in the
    % leaf model
end
% [unfinished] check prep_bnd(:,3). Interior boundaries should not have -1
prep_bnd(find(prep_bnd(:,3)==-1),3)=0;%currently just let those=0, i.e. set flux=0 on those boundaries

%% create selection manually
model.selection.create('bnd', 'Explicit');
model.selection('bnd').geom('geom1', 2);
model.selection('bnd').set(prep_bnd(find(prep_bnd(:,3)==0),1)');
model.selection('bnd').name('model boundary');
model.selection.create('air', 'Explicit');
model.selection('air').geom('geom1', 3);
model.selection('air').set(prep_bnd(find(prep_dom(:,3)==5),1)');
model.selection('air').name('air');
model.selection.create('mscyto', 'Explicit');
model.selection('mscyto').geom('geom1', 3);
model.selection('mscyto').set(prep_bnd(find(prep_dom(:,3)==1),1)');
model.selection('mscyto').name('ms cyto');
model.selection.create('mschlo', 'Explicit');
model.selection('mschlo').geom('geom1', 3);
model.selection('mschlo').set(prep_bnd(find(prep_dom(:,3)==2),1)');
model.selection('mschlo').name('ms chlo');
model.selection.create('msmito', 'Explicit');
model.selection('msmito').geom('geom1', 3);
model.selection('msmito').set(prep_bnd(find(prep_dom(:,3)==3),1)');
model.selection('msmito').name('ms mito');
model.selection.create('msvacu', 'Explicit');
model.selection('msvacu').geom('geom1', 3);
model.selection('msvacu').set(prep_bnd(find(prep_dom(:,3)==4),1)');
model.selection('msvacu').name('ms vacu');
model.selection.create('bscyto', 'Explicit');
model.selection('bscyto').geom('geom1', 3);
model.selection('bscyto').set(prep_bnd(find(prep_dom(:,3)==8),1)');
model.selection('bscyto').name('bs cyto');
model.selection.create('bschlo', 'Explicit');
model.selection('bschlo').geom('geom1', 3);
model.selection('bschlo').set(prep_bnd(find(prep_dom(:,3)==9),1)');
model.selection('bschlo').name('bs chlo');
model.selection.create('bsmito', 'Explicit');
model.selection('bsmito').geom('geom1', 3);
model.selection('bsmito').set(prep_bnd(find(prep_dom(:,3)==10),1)');
model.selection('bsmito').name('bs mito');
model.selection.create('stom', 'Explicit');
model.selection('stom').geom('geom1', 3);
model.selection('stom').set(prep_bnd(find(prep_dom(:,3)==6|prep_dom(:,3)==7),1)');
model.selection('stom').name('stom');
model.selection.create('uni1', 'Union');
model.selection('uni1').set('input', {'bschlo' 'mschlo'});
model.selection('uni1').label('all chlo');


%%volume of all chloroplast
tmp_sel=mphgetselection(model.selection('bschlo'));
tmp_bnd_num=tmp_sel.entities;
model.geom('geom1').measureFinal.selection.geom('geom1', 3);
model.geom('geom1').measureFinal.selection.set(tmp_bnd_num);
bsc_chlo_vol=model.geom('geom1').measureFinal().getVolume();
tmp_sel=mphgetselection(model.selection('mschlo'));
tmp_bnd_num=tmp_sel.entities;
model.geom('geom1').measureFinal.selection.geom('geom1', 3);
model.geom('geom1').measureFinal.selection.set(tmp_bnd_num);
msc_chlo_vol=model.geom('geom1').measureFinal().getVolume();
all_chlo_vol=bsc_chlo_vol+msc_chlo_vol;
xmax=BSE_width/2+EPU_width+BUT_width/2;
ymax=(EPL_thick+MST_thickatvein+EPU_thick)+0.1e-6;
zmax=LEAF_z;
xmin=0;ymin=-0.1e-6;zmin=0;
leaf_surface_area=xmax*zmax;
AVR=all_chlo_vol/leaf_surface_area/2;

%convert [chl] from mg mm-2 to g/m3
% ([chl]*1e-3*1e6)/(AVR*2)
%chl_con=0.000664;%mg/mm2
%chl_con4raytracing=chl_con*1e3/(AVR*2)

save save_e_geom.mat -regexp '^(?!(model|ans)$).'

