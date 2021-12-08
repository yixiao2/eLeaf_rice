% eLeaf: 3D model of rice leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 1.2.6

tic
%path('C:\Program Files\COMSOL\COMSOL52a\Multiphysics\mli',path)
%mphstart(2036);
import com.comsol.model.*
import com.comsol.model.util.*
model=mphload('tmpCK_geomIP_cad_mesh_cresel.mph');
load save_e_geom.mat %from e_geo_main_v1_2.m
% copyfile('../2.e_raytracing/results_cell','results_cell');

e_pre_physics_v1_2_3
load ab_profile.mat

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

%% create integration
model.cpl.create('intop1', 'Integration', 'geom1');
model.cpl('intop1').selection.named('mschlo');
model.cpl('intop1').label('mschlo');
model.cpl.create('intop2', 'Integration', 'geom1');
model.cpl('intop2').selection.named('bschlo');
model.cpl('intop2').label('bschlo');
model.cpl.create('intop3', 'Integration', 'geom1');
model.cpl('intop3').selection.named('msmito');
model.cpl('intop3').label('msmito');
model.cpl.create('intop4', 'Integration', 'geom1');
model.cpl('intop4').selection.named('bsmito');
model.cpl('intop4').label('bsmito');
for tmp_loop=1:count_chl
    tmp_str1=['mschlo_',num2str(tmp_loop)];
    tmp_str2=['intop',num2str(4+tmp_loop*2-1)];
    model.cpl.create(tmp_str2, 'Integration', 'geom1');
    model.cpl(tmp_str2).selection.named(tmp_str1);
    model.cpl(tmp_str2).label(tmp_str1);
    tmp_str1=['msmito_',num2str(tmp_loop)];
    tmp_str2=['intop',num2str(4+tmp_loop*2)];
    model.cpl.create(tmp_str2, 'Integration', 'geom1');
    model.cpl(tmp_str2).selection.named(tmp_str1);
    model.cpl(tmp_str2).label(tmp_str1);
end

%% global variable and model variable
count_var=0;
%%%bs
% bs cyto
model.variable.create('var1');count_var=count_var+1;
model.variable('var1').model('mod1');
model.variable('var1').set('visc', 'Vcy');%visc is short for viscosity
model.variable('var1').set('Xa', 'Xac');
model.variable('var1').set('pH', 'pHc');
model.variable('var1').set('hy', 'ka*Xa*(c1-c2/(Keq/(10^-pH)))/(KmCA+KmCA/KmBA*c2+c1)');
model.variable('var1').set('fvc', '0');%fvc is the carboxylation flux rate
model.variable('var1').set('fphres', '0');%fph is the photorespiration flux rate
model.variable('var1').set('fres', '0');%fres is the respiration rate flux rate
model.variable('var1').set('c1_ini','ppi*sc');
model.variable('var1').set('c2_ini','ppi*sc*10');
model.variable('var1').selection.geom('geom1', 3);
model.variable('var1').selection.named('bscyto');
model.variable('var1').label('bs cyto');
% bs chlo
model.variable.create('var2');count_var=count_var+1;
model.variable('var2').model('mod1');
model.variable('var2').set('visc', 'Vchl');%visc is short for viscosity
model.variable('var2').set('Xa', 'Xas');
model.variable('var2').set('pH', 'pHs');
model.variable('var2').set('hy', 'ka*Xa*(c1-c2/(Keq/(10^-pH)))/(KmCA+KmCA/KmBA*c2+c1)');
model.variable('var2').set('fc', 'Vcmax_fvcb*leafsurface/Vols_all*c1/(c1+Km_fvcb)', 'CO2 lim');
model.variable('var2').set('fj', 'J*c1/(4*c1+8*K_fvcb)', 'RuBP lim');
model.variable('var2').set('fvc', 'min(fc,fj)');%fvc is the carboxylation flux rate
model.variable('var2').set('fvc_bk', 'co2fix_eps_83m100(c1,J_eps)');%v1.2
model.variable('var2').set('fphres', '0');%fph is the photorespiration flux rate
model.variable('var2').set('fres', '0');%fres is the respiration rate flux rate
model.variable('var2').set('integ_src_fphres', 'fvc/c1*K_fvcb');
model.variable('var2').set('integ_src_fphres_bk', 'co2rsp_eps_83m100(c1,J_eps)');%v1.2
model.variable('var2').set('Vols_bs','mod1.intop2(1)');
model.variable('var2').set('Jmax', 'Jmax_fvcb*leafsurface/Vols_all*Vols_bs');
tmp_ab=ab_profile(1);
model.variable('var2').set('ab_bs',num2str(tmp_ab));
model.variable('var2').set('I2', 'irra*leafsurface*ab_bs*(s_yin*yiill/ab_all)');
model.variable('var2').set('J', '(I2+Jmax-sqrt((I2+Jmax)^2-4*theta_fvcb*I2*Jmax))/(2*theta_fvcb)/Vols_bs');
model.variable('var2').set('J_eps','irra*leafsurface*ab_bs/Vols_bs');
model.variable('var2').set('c1_ini','ppi*sc*0.7');
model.variable('var2').set('c2_ini','ppi*sc*40');
model.variable('var2').selection.geom('geom1', 3);
model.variable('var2').selection.named('bschlo');
model.variable('var2').label('bs chlo');
% bs mito
model.variable.create('var3');count_var=count_var+1;
model.variable('var3').model('mod1');
model.variable('var3').set('visc', 'Vcy');%visc is short for viscosity
model.variable('var3').set('Xa', 'Xac');
model.variable('var3').set('pH', 'pHc');
model.variable('var3').set('hy', 'ka*Xa*(c1-c2/(Keq/(10^-pH)))/(KmCA+KmCA/KmBA*c2+c1)');
model.variable('var3').set('fvc', '0');%fvc is the carboxylation flux rate
model.variable('var3').set('Volm_bs','mod1.intop4(1)');
model.variable('var3').set('fphres', 'integ_fphres/Volm_bs');%fph is the photorespiration flux rate
model.variable('var3').set('fres', 'Resp*leafsurface/Volm_all');%fres is the respiration rate flux rate
model.variable('var3').set('integ_fphres', 'mod1.intop2(mod1.integ_src_fphres)', 'total amount of oxygenation');
model.variable('var3').set('c1_ini','ppi*sc*1.5');
model.variable('var3').set('c2_ini','ppi*sc*40');
model.variable('var3').selection.geom('geom1', 3);
model.variable('var3').selection.named('bsmito');
model.variable('var3').label('bs mito');
%%%ms
% ms cyto
model.variable.create('var4');count_var=count_var+1;
model.variable('var4').model('mod1');
model.variable('var4').set('visc', 'Vcy');%visc is short for viscosity
model.variable('var4').set('Xa', 'Xac');
model.variable('var4').set('pH', 'pHc');
model.variable('var4').set('hy', 'ka*Xa*(c1-c2/(Keq/(10^-pH)))/(KmCA+KmCA/KmBA*c2+c1)');
model.variable('var4').set('fvc', '0');%fvc is the carboxylation flux rate
model.variable('var4').set('fphres', '0');%fph is the photorespiration flux rate
model.variable('var4').set('fres', '0');%fres is the respiration rate flux rate
model.variable('var4').set('c1_ini','ppi*sc');
model.variable('var4').set('c2_ini','ppi*sc*10');
model.variable('var4').selection.geom('geom1', 3);
model.variable('var4').selection.named('mscyto');
model.variable('var4').label('ms cyto');
% ms vacu
model.variable.create('var5');count_var=count_var+1;
model.variable('var5').model('mod1');
model.variable('var5').set('visc', 'Vcy');%visc is short for viscosity
%model.variable('var5').set('Xa', 'Xac');
%model.variable('var5').set('pH', 'pHc');
model.variable('var5').set('hy', '0');
model.variable('var5').set('fvc', '0');%fvc is the carboxylation flux rate
model.variable('var5').set('fphres', '0');%fph is the photorespiration flux rate
model.variable('var5').set('fres', '0');%fres is the respiration rate flux rate
model.variable('var5').set('c1_ini','ppi*sc');
model.variable('var5').set('c2_ini','ppi*sc*10');
model.variable('var5').selection.geom('geom1', 3);
model.variable('var5').selection.named('msvacu');
model.variable('var5').label('ms vacu');
% chlo
for tmp_loop=1:count_chl
    count_var=count_var+1;tmp_str=['var',num2str(count_var)];model.variable.create(tmp_str);
    tmp_str1=['mschlo_',num2str(tmp_loop)];
    model.variable(tmp_str).model('mod1');
    model.variable(tmp_str).set('visc', 'Vchl');%visc is short for viscosity
    model.variable(tmp_str).set('Xa', 'Xas');
    model.variable(tmp_str).set('pH', 'pHs');
    model.variable(tmp_str).set('hy', 'ka*Xa*(c1-c2/(Keq/(10^-pH)))/(KmCA+KmCA/KmBA*c2+c1)');
    model.variable(tmp_str).set('fc', 'Vcmax_fvcb*leafsurface/Vols_all*c1/(c1+Km_fvcb)', 'CO2 lim');
    model.variable(tmp_str).set('fj', 'J*c1/(4*c1+8*K_fvcb)', 'RuBP lim');
    model.variable(tmp_str).set('fvc', 'min(fc,fj)');%fvc is the carboxylation flux rate
    model.variable(tmp_str).set('fvc_bk', 'co2fix_eps_83m100(c1,J_eps)');%v1.2
    model.variable(tmp_str).set('fphres', '0');%fph is the photorespiration flux rate
    model.variable(tmp_str).set('fres', '0');%fres is the respiration rate flux rate
    model.variable(tmp_str).set('integ_src_fphres', 'fvc/c1*K_fvcb');
    model.variable(tmp_str).set('integ_src_fphres_bk', 'co2rsp_eps_83m100(c1,J_eps)');%v1.2
    %model.variable(tmp_str).set('Jmax', 'Jmax_leaf*leafsurface/Vols*Vols_bs');
    %model.variable(tmp_str).set('I2', 'irra*leafsurface*2*ab_bs*(1-0.15)/2');
    %model.variable(tmp_str).set('J', '(I2+Jmax-sqrt((I2+Jmax)^2-4*theta*I2*Jmax))/(2*theta)/Vols_bs');
    tmp_str2=['mod1.intop',num2str(4+tmp_loop*2-1),'(1)'];%see above line165>>tmp_str2=['intop',num2str(4+tmp_loop*2-1)];
    model.variable(tmp_str).set('Vols_1ms',tmp_str2);
    model.variable(tmp_str).set('Jmax', 'Jmax_fvcb*leafsurface/Vols_all*Vols_1ms');
    tmp_ab=ab_profile(tmp_loop+1);%**link with ray tracing results**%
    model.variable(tmp_str).set('ab_1ms',num2str(tmp_ab));
    model.variable(tmp_str).set('I2', 'irra*leafsurface*ab_1ms*(s_yin*yiill/ab_all)');
    model.variable(tmp_str).set('J', '(I2+Jmax-sqrt((I2+Jmax)^2-4*theta_fvcb*I2*Jmax))/(2*theta_fvcb)/Vols_1ms');
    model.variable(tmp_str).set('J_eps','irra*leafsurface*ab_1ms/Vols_1ms');
    model.variable(tmp_str).set('c1_ini','ppi*sc*0.7');
    model.variable(tmp_str).set('c2_ini','ppi*sc*40');
    model.variable(tmp_str).selection.geom('geom1', 3);
    model.variable(tmp_str).selection.named(tmp_str1);
    model.variable(tmp_str).label(tmp_str1);
end
% mito
for tmp_loop=1:count_chl
    count_var=count_var+1;tmp_str=['var',num2str(count_var)];model.variable.create(tmp_str);
    tmp_str1=['msmito_',num2str(tmp_loop)];
    model.variable(tmp_str).model('mod1');
    model.variable(tmp_str).set('visc', 'Vcy');%visc is short for viscosity
    model.variable(tmp_str).set('Xa', 'Xac');
    model.variable(tmp_str).set('pH', 'pHc');
    model.variable(tmp_str).set('hy', 'ka*Xa*(c1-c2/(Keq/(10^-pH)))/(KmCA+KmCA/KmBA*c2+c1)');
    model.variable(tmp_str).set('fvc', '0');%fvc is the carboxylation flux rate
%     model.variable(tmp_str).set('fphres', 'integ_fphres/Volm_bs');%fph is the photorespiration flux rate
%     model.variable(tmp_str).set('fres', 'Resp*leafsurface/Volm_all');%fres is the respiration rate flux rate
%     model.variable(tmp_str).set('integ_fphres', 'mod1.intop2(mod1.integ_src_fphres)', 'total amount of oxygenation');
    tmp_str2=['mod1.intop',num2str(4+tmp_loop*2),'(1)'];%see line170
    model.variable(tmp_str).set('Volm_1ms',tmp_str2);
    model.variable(tmp_str).set('fphres', 'integ_fphres/Volm_1ms');
    model.variable(tmp_str).set('fres', 'Resp*leafsurface/Volm_all');
    tmp_str3=['mod1.intop',num2str(4+tmp_loop*2-1),'(mod1.integ_src_fphres)'];
    model.variable(tmp_str).set('integ_fphres',tmp_str3, 'total amount of oxygenation');
    model.variable(tmp_str).set('c1_ini','ppi*sc*1.5');
    model.variable(tmp_str).set('c2_ini','ppi*sc*40');
    model.variable(tmp_str).selection.geom('geom1', 3);
    model.variable(tmp_str).selection.named(tmp_str1);
    model.variable(tmp_str).label(tmp_str1);
end
% air
count_var=count_var+1;tmp_str=['var',num2str(count_var)];model.variable.create(tmp_str);
model.variable(tmp_str).model('mod1');
model.variable(tmp_str).set('visc', 'Vair');
model.variable(tmp_str).set('hy', '0');
model.variable(tmp_str).set('fvc', '0');%fvc is the carboxylation flux rate
model.variable(tmp_str).set('fphres', '0');%fph is the photorespiration flux rate
model.variable(tmp_str).set('fres', '0');%fres is the respiration rate flux rate
model.variable(tmp_str).set('c1_ini','ppi*sc');
model.variable(tmp_str).set('c2_ini','0');
model.variable(tmp_str).selection.geom('geom1', 3);
model.variable(tmp_str).selection.named('air');
model.variable(tmp_str).label('air');
% stomata
count_var=count_var+1;tmp_str=['var',num2str(count_var)];model.variable.create(tmp_str);
model.variable(tmp_str).model('mod1');
model.variable(tmp_str).set('visc', 'Vair');
model.variable(tmp_str).set('hy', '0');
model.variable(tmp_str).set('fvc', '0');%fvc is the carboxylation flux rate
model.variable(tmp_str).set('fphres', '0');%fph is the photorespiration flux rate
model.variable(tmp_str).set('fres', '0');%fres is the respiration rate flux rate
model.variable(tmp_str).set('c1_ini','ppi*sc');
model.variable(tmp_str).set('c2_ini','0');
model.variable(tmp_str).selection.geom('geom1', 3);
model.variable(tmp_str).selection.named('stom');
model.variable(tmp_str).label('stom');

% model variable on entire domain
count_var=count_var+1;tmp_str=['var',num2str(count_var)];model.variable.create(tmp_str);
model.variable(tmp_str).set('Vols_all','mod1.intop1(1)+mod1.intop2(1)');
model.variable(tmp_str).set('Volm_all','mod1.intop3(1)+mod1.intop4(1)');
% model.variable(tmp_str).set('Vols','mod1.intop2(1)');
% model.variable(tmp_str).set('Volm','mod1.intop3(1)');
% model.variable(tmp_str).set('Volv','mod1.intop4(1)');
% model.variable(tmp_str).set('Volc','mod1.intop1(1)');

% additional model parameters
model.param.set('Dc', '1.83e-9[m^2/s]', 'diff coef co2');
model.param.set('Db', '0.52*Dc', 'dif cons HCO3-');
model.param.set('R', '8.31[J/K/mol]', 'gas constant');
%model.param.set('Xa', '0.3[mol/m^3]', 'CA concentration');
model.param.set('Xa', [num2str(0.3*fc_CA),'[mol/m^3]'], 'CA concentration');
%model.param.set('Xa2', '0.5[mol/m^3]', 'CA concentration cytosol');
model.param.set('Xa2', [num2str(0.15*fc_CA),'[mol/m^3]'], 'CA concentration');
model.param.set('Xac', 'Xa2', 'CA concentration cytosol');
model.param.set('Xas', 'Xa', 'CA concentration cytosol');
model.param.set('Xam', 'Xa', 'CA concentration cytosol');
model.param.set('ka', '300000[1/s]', 'kcat carbonic anhydrase 25degC');
model.param.set('KmCA', '1.5[mol/m^3]', 'Km CA at 25degC');
%model.param.set('dMSwall','0.16[um]');
model.param.set('dMSwall',[num2str(cellwallthick),'[um]']);
model.param.set('DcMS','Dc*ptMS');
%model.param.set('ptMS','0.2','effective porosity MS wall');
model.param.set('ptMS',num2str(0.2*fc_wall_perm),'effective porosity MS wall');
%model.param.set('dBSwall','0.16[um]');
model.param.set('dBSwall',[num2str(cellwallthick),'[um]']);
model.param.set('DcBS','Dc*ptBS');
%model.param.set('ptBS','0.1','effective porosity BS wall');
model.param.set('ptBS',num2str(0.1*fc_wall_perm),'effective porosity BS wall');
model.param.set('pHs', '8.0', 'stroma pH');
model.param.set('pO2', '21[kPa]', 'oxygen concentration');
%model.param.set('Resp', '0.283[umol/m^2/s]', 'dark respiration in light');
model.param.set('Resp', [num2str(Resp),'[umol/m^2/s]'],'dark respiration in light');
model.param.set('ppi_bak', '30[Pa]', 'internal co2 conc');
model.param.set('Dmenv', 'Pco2_m*tmenv', 'mito enveloppe diff const');
model.param.set('tmenv', '0.02e-6[m]', 'mito enveloppe thickness');
model.param.set('Dcenv', 'Pco2*tcenv', 'chloro envenloppe diff constant');
model.param.set('tcenv', 'tmenv', 'chloro enveloppe thickness');
model.param.set('stdT','273.15[K]+25[K]');
model.param.set('R','8.31[J/K/mol]');
model.param.set('Dair_true','0.1381[cm^2/s]*(stdT/282[K])^(1.81)');%% Massman 1998
model.param.set('Vair','Dc/Dair*sc*R*stdT');
model.param.set('Vchl', '10', 'relative viscosity chloroplasts');
model.param.set('Vmit', 'Vchl', 'relative viscos mitochondria');
model.param.set('Vcy', '2', 'relative viscosity cytosol');
model.param.set('pHm', '8', 'pH mitochondria');
model.param.set('Dcenvb', 'Phco3*tcenv', 'hco3 diffusion constant chl mem');
model.param.set('KmBA', '34[mol/m^3]', 'Michaelis constant CA');
model.param.set('Keq', '5.6e-7', 'equilibrium constant CA');
model.param.set('Pco2', '0.35[cm/s]', 'chl mem permeability co2');
model.param.set('Pco2_m', '0.35[cm/s]', 'mit mem permeability co2');
model.param.set('Phco3', '5e-7[m/s]', 'chl mem permeability hco3');
model.param.set('Phco3_m', '5e-7[m/s]', 'mit mem permeability hco3');
model.param.set('Dmenvb', 'Phco3_m*tmenv', 'hco3 diffusion constant mit mem');
model.param.set('K_bak', '1.35e-3[mol/m^3]', 'Gamma-star');
model.param.set('sc', '0.33e-3[mol/m^3/Pa]', 'solubility co2');%% c=p*sc
model.param.set('Pair', '100[kPa]', 'air pressure');
model.param.set('pHc', '7.3', 'pH cytosol');
%model.param.set('leafsurface2', '(1.22e-5*2*5e-6*2)[m^2]', 'my estimate of model leaf area');
model.param.set('leafsurface', '(BSE_width/2+EPU_width+BUT_width/2)*LEAF_z');
model.param.set('v1_stable', '1.2336648[mol/(m^3*s)]');
model.param.set('so', '1.26[mol/m^3/bar]', 'oxygen solubility');
model.param.set('Ci_gas', '0.01[mol/m^3]');
model.param.set('Ci','Ci_gas*R*stdT*sc');%%trick to handle gas phase and liquid phase;
model.param.set('Dair','Dair_true/R/stdT/sc');%% Massman 1998
%% [Causion], output C --> convert to partial pressure, C/sc. Partial pressure is continous. So concentration of true gaseous Ci_gas = Ci/sc/RT
model.param.set('ppi', 'Ci/sc');
model.param.set('irra', '1000[umol/m^2/s]');
model.param.set('O2', 'pO2*so');
%model.param.set('Vcmax_fvcb', '191.7[umol/m^2/s]');
model.param.set('Vcmax_fvcb', [num2str(Vcmax),'[umol/m^2/s]']);
%model.param.set('Sc_o_fvcb', '2138', 'Sc/o (bar/bar)');
model.param.set('Sc_o_fvcb', num2str(Sco), 'Sc/o (bar/bar)');
model.param.set('K_pp_fvcb', '0.5*pO2/Sc_o_fvcb');
model.param.set('K_fvcb', 'K_pp_fvcb*sc');
%model.param.set('Kc_fvcb', '23.9[Pa]*sc');
model.param.set('Kc_fvcb', [num2str(Kc/10),'[Pa]*sc']);
%model.param.set('Ko_fvcb', '266[mbar]*so');
model.param.set('Ko_fvcb', [num2str(Ko/1000),'[mbar]*so']);
model.param.set('Km_fvcb', 'Kc_fvcb*(1+O2/Ko_fvcb)');
%model.param.set('Jmax_fvcb', '249.7[umol/m^2/s]');
model.param.set('Jmax_fvcb', [num2str(Jm),'[umol/m^2/s]']);
%model.param.set('theta_fvcb', '0.98');
model.param.set('theta_fvcb', num2str(theta));
%model.param.set('s_yin', '0.572');
model.param.set('s_yin', num2str(s_yin));
%model.param.set('yiill','0.653');
model.param.set('yiill',num2str(yiill));
model.param.set('ab_all',num2str(ab_all));


%% physics
model.physics.create('chds', 'DilutedSpecies', 'geom1');
%model.physics('chds').field('concentration').field('c1');
model.physics('chds').field('concentration').component({'c1' 'c2'});
%model.physics('chds').prop('Convection').set('Convection', '0');
model.physics('chds').prop('TransportMechanism').set('Convection', false);
model.physics('chds').feature('cdm1').set('D_c1', {'Dc/visc' '0' '0' '0' 'Dc/visc' '0' '0' '0' 'Dc/visc'});
model.physics('chds').feature('cdm1').set('D_c2', {'Db/visc' '0' '0' '0' 'Db/visc' '0' '0' '0' 'Db/visc'});

model.physics('chds').feature.create('reac1', 'Reactions', 3);
model.physics('chds').feature('reac1').selection.all;
%model.physics('chds').feature('reac1').set('R', {'-hy'; 'hy'});
model.physics('chds').feature('reac1').setIndex('R_c1', '-(hy+fvc)+fphres+fres', 0);
model.physics('chds').feature('reac1').setIndex('R_c2', 'hy', 0);

model.physics('chds').feature('init1').setIndex('initc', 'c1_ini', 0);
model.physics('chds').feature('init1').setIndex('initc', 'c2_ini', 1);

%for stomata Ci
model.physics('chds').create('conc1', 'PairConcentration', 2);
tmp_pairString=TagsPairs(prep_pairs(find(prep_pairs(:,2)==6|prep_pairs(:,2)==7),1));
clear tmp_String;
tmp_String=('''');
for tmp_loop=1:size(tmp_pairString,1)
    tmp_String=[tmp_String,char(tmp_pairString(tmp_loop,1)),''''];
    if tmp_loop~=size(tmp_pairString,1)
        tmp_String=[tmp_String,';',''''];
    end
end
eval(['model.physics(''chds'').feature(''conc1'').set(''pairs'', {',tmp_String,'});'])
model.physics('chds').feature('conc1').setIndex('species', true, 0);
model.physics('chds').feature('conc1').setIndex('c0', 'Ci', 0);
% %for IAS
% %model.physics('chds').feature.create('cdm2', 'ConvectionDiffusion',3);%v43b
% model.physics('chds').create('cdm2', 'ConvectionDiffusionMigration', 3);%v52a
% model.physics('chds').feature('cdm2').selection.named('air');
% model.physics('chds').feature('cdm2').set('D_c1', {'Dc/visc'; '0'; '0'; '0'; 'Dc/visc'; '0'; '0'; '0'; 'Dc/visc'});
% model.physics('chds').feature('cdm2').set('D_c2', {'Db/visc'; '0'; '0'; '0'; 'Db/visc'; '0'; '0'; '0'; 'Db/visc'});
% model.physics('chds').feature('cdm2').name('Diffusion 2');
%for IAS-MS
model.physics('chds').feature.create('ptdb1', 'PairThinDiffusionBarrier', 2);
tmp_pairString=TagsPairs(prep_pairs(find(prep_pairs(:,2)==5),1));
clear tmp_String;
tmp_String=('''');
for tmp_loop=1:size(tmp_pairString,1)
    tmp_String=[tmp_String,char(tmp_pairString(tmp_loop,1)),''''];
    if tmp_loop~=size(tmp_pairString,1)
        tmp_String=[tmp_String,';',''''];
    end
end
eval(['model.physics(''chds'').feature(''ptdb1'').set(''pairs'', {',tmp_String,'});'])
model.physics('chds').feature('ptdb1').set('Ds', {'DcMS'; '0'});
model.physics('chds').feature('ptdb1').set('ds', 'dMSwall');
%for IAS-BS
model.physics('chds').feature.create('ptdb2', 'PairThinDiffusionBarrier', 2);
tmp_pairString=TagsPairs(prep_pairs(find(prep_pairs(:,2)==8),1));
clear tmp_String;
tmp_String=('''');
for tmp_loop=1:size(tmp_pairString,1)
    tmp_String=[tmp_String,char(tmp_pairString(tmp_loop,1)),''''];
    if tmp_loop~=size(tmp_pairString,1)
        tmp_String=[tmp_String,';',''''];
    end
end
eval(['model.physics(''chds'').feature(''ptdb2'').set(''pairs'', {',tmp_String,'});'])
model.physics('chds').feature('ptdb2').set('Ds', {'DcBS'; '0'});
model.physics('chds').feature('ptdb2').set('ds', 'dBSwall');
%for MS-BS
model.physics('chds').feature.create('ptdb3', 'PairThinDiffusionBarrier', 2);
tmp_pairString=TagsPairs(prep_pairs(find(prep_pairs(:,2)==9),1));
clear tmp_String;
tmp_String=('''');
for tmp_loop=1:size(tmp_pairString,1)
    tmp_String=[tmp_String,char(tmp_pairString(tmp_loop,1)),''''];
    if tmp_loop~=size(tmp_pairString,1)
        tmp_String=[tmp_String,';',''''];
    end
end
if isempty(tmp_pairString)~=1
    eval(['model.physics(''chds'').feature(''ptdb3'').set(''pairs'', {',tmp_String,'});'])
    model.physics('chds').feature('ptdb3').set('Ds', {'1/((dMSwall/(dMSwall+dBSwall)/DcMS)+(dBSwall/(dMSwall+dBSwall)/DcBS))'; '0'});
    model.physics('chds').feature('ptdb3').set('ds', 'dMSwall+dBSwall');
end
%for MS-MS
model.physics('chds').feature.create('ptdb4', 'PairThinDiffusionBarrier', 2);
tmp_pairString=TagsPairs(prep_pairs(find(prep_pairs(:,2)==1),1));
clear tmp_String;
tmp_String=('''');
for tmp_loop=1:size(tmp_pairString,1)
    tmp_String=[tmp_String,char(tmp_pairString(tmp_loop,1)),''''];
    if tmp_loop~=size(tmp_pairString,1)
        tmp_String=[tmp_String,';',''''];
    end
end
eval(['model.physics(''chds'').feature(''ptdb4'').set(''pairs'', {',tmp_String,'});'])
model.physics('chds').feature('ptdb4').set('Ds', {'DcMS'; '0'});
model.physics('chds').feature('ptdb4').set('ds', 'dMSwall*2');
%for cyto-chlo
model.physics('chds').feature.create('ptdb5', 'PairThinDiffusionBarrier', 2);
tmp_pairString=TagsPairs(prep_pairs(find(prep_pairs(:,2)==2|prep_pairs(:,2)==10),1));
clear tmp_String;
tmp_String=('''');
for tmp_loop=1:size(tmp_pairString,1)
    tmp_String=[tmp_String,char(tmp_pairString(tmp_loop,1)),''''];
    if tmp_loop~=size(tmp_pairString,1)
        tmp_String=[tmp_String,';',''''];
    end
end
eval(['model.physics(''chds'').feature(''ptdb5'').set(''pairs'', {',tmp_String,'});'])
model.physics('chds').feature('ptdb5').set('Ds', {'Dcenv'; 'Dcenvb'});
model.physics('chds').feature('ptdb5').set('ds', 'tcenv');
%for cyto-mito
model.physics('chds').feature.create('ptdb6', 'PairThinDiffusionBarrier', 2);
tmp_pairString=TagsPairs(prep_pairs(find(prep_pairs(:,2)==3|prep_pairs(:,2)==11),1));
clear tmp_String;
tmp_String=('''');
for tmp_loop=1:size(tmp_pairString,1)
    tmp_String=[tmp_String,char(tmp_pairString(tmp_loop,1)),''''];
    if tmp_loop~=size(tmp_pairString,1)
        tmp_String=[tmp_String,';',''''];
    end
end
eval(['model.physics(''chds'').feature(''ptdb6'').set(''pairs'', {',tmp_String,'});'])
model.physics('chds').feature('ptdb6').set('Ds', {'Dmenv'; 'Dmenvb'});
model.physics('chds').feature('ptdb6').set('ds', 'tmenv');
%for cyto-vacu
model.physics('chds').feature.create('ptdb7', 'PairThinDiffusionBarrier', 2);
tmp_pairString=TagsPairs(prep_pairs(find(prep_pairs(:,2)==4),1));
clear tmp_String;
tmp_String=('''');
for tmp_loop=1:size(tmp_pairString,1)
    tmp_String=[tmp_String,char(tmp_pairString(tmp_loop,1)),''''];
    if tmp_loop~=size(tmp_pairString,1)
        tmp_String=[tmp_String,';',''''];
    end
end
eval(['model.physics(''chds'').feature(''ptdb7'').set(''pairs'', {',tmp_String,'});'])
model.physics('chds').feature('ptdb7').set('Ds', {'Dmenv'; 'Dmenvb'});
model.physics('chds').feature('ptdb7').set('ds', 'tmenv');

% model.mesh.create('mesh1', 'geom1');
% model.mesh('mesh1').feature.create('ftet1', 'FreeTet');
% model.mesh('mesh1').run;

if ARG_TYPE==0
    mphsave(model,'eleaf_fvcb_CKIR64.mph');
else
    mphsave(model,'eleaf_fvcb_HCIR64.mph');
end

toc
%end
