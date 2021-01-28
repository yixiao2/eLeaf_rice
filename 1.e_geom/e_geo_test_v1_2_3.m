% eLeaf: 3D model of rice leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 1.2.5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update 2017/12/27
% - mph geometry pass the test of free tetrahedron mesh
% before export object for ray tracing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model.geom('geom1').feature('fin').name('Form Assembly');
model.geom('geom1').feature('fin').set('action', 'assembly');
model.geom('geom1').feature('fin').set('createpairs', 'on');
model.geom('geom1').feature('fin').set('imprint', 'on');
%model.geom('geom1').feature('fin').set('repairtol', '1e-5');
%model.geom('geom1').feature(tmptag_IAS).set('repairtoltype', 'auto');
model.geom('geom1').run;
mphsave(model,'tmpCK_geomIP_cad_nocresel.mph')

%% test free tetrahedron
try
    model.mesh.create('mesh1', 'geom1');
    model.mesh('mesh1').feature.create('ftet1', 'FreeTet');
    model.mesh('mesh1').run;
    
    % output tags for create selection
    if loop_idx==-1&&loop_idy==-1
        save save_e_geom.mat -regexp '^(?!(model|ans)$).'
    end
    save list_cresel.mat list_cresel
    toc
catch
    display('FAIL free tetrahedron mesh');
    ModelUtil.remove('Model_IAS2D');
    ModelUtil.remove('Model_3Dleaf')
    e_geo_IAS2D_v1_2_1
    e_geo_main_v1_2_1
    e_geo_test_v1_2_3
end
mphsave(model,'tmpCK_geomIP_cad_mesh_nocresel.mph')
%% end test

