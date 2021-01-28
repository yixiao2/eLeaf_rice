% eLeaf: 3D model of rice leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 1.2.5

%% step1: get eleaf_fvcb_{CKIR64,HCIR64}_a_phipsii.mph
%ARG_TYPE=input('ARG_TYPE (CK=0; HC=1):');
e_physics_v1_2_3_fvcb;%eleaf_fvcb_{CKIR64,HCIR64}.mph
e_study_v1_2_3_prswp_forfvcb_fit;%eleaf_fvcb_{CKIR64,HCIR64}_a_phipsii.mph
