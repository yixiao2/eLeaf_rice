% eLeaf: 3D model of rice leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 1.2.5

function e_physics_addcresel_v1_1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update 2017/12/27
% - mph geometry pass the test of free tetrahedron mesh
% before export object for ray tracing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%path('C:\Program Files\COMSOL\COMSOL52a\Multiphysics\mli',path)
%mphstart(2036);
%import com.comsol.model.*
%import com.comsol.model.util.*
model=mphload('tmpCK_geomIP_cad_mesh_nocresel.mph');
load list_cresel
for count_list_cresel=1:size(list_cresel,2)
    tmptag_cresel=list_cresel(count_list_cresel);
    model.geom('geom1').feature(tmptag_cresel).set('createselection', true);
end
model.geom('geom1').run;
mphsave(model,'tmpCK_geomIP_cad_mesh_cresel.mph');
%ModelUtil.disconnect;
%exit;
end
