% eLeaf: 3D model of rice leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <xiaoyi@sippe.ac.cn>
% @version: 1.2.5

%ARG_TYPE=input('ARG_TYPE/model name (CK=0; HC=1):');
function run_e_leaf_fvcb_SA()
ARG_TYPE=0;
e_physics_v1_2_3_fvcb;
e_study_v1_2_5_prswp_forfvcb_dissection;
end
