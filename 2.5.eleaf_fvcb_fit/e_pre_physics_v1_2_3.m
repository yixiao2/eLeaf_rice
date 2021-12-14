% eLeaf: 3D model of rice leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 1.2.5

load save_e_geom.mat %from e_geo_main_v1_2.m
copyfile('../2.e_raytracing/results_merged*','./');

%% readin results from ray tracing: 475nm and 625nm
load save_raytracing.mat count_mschl
tmp_ab_profile=zeros(2,size(count_mschl,2));

tmp_rtsum=importdata('results_merged_rtsum_475nm_500x_rep1');
tmp_ab_chlo_results=importdata('results_merged_absrf_475nm_500x_rep1');
tmp_ab_profile(1,:)=tmp_ab_chlo_results(tmp_ab_chlo_results(:,3)==1,7)/(tmp_rtsum.data(5));
spectrum_retrab(1,1)=tmp_rtsum.data(2)/tmp_rtsum.data(5);%re;
spectrum_retrab(1,2)=tmp_rtsum.data(3)/tmp_rtsum.data(5);%re;
spectrum_retrab(1,3)=tmp_rtsum.data(1)/tmp_rtsum.data(5);%re;

tmp_rtsum=importdata('results_merged_rtsum_625nm_500x_rep1');
tmp_ab_chlo_results=importdata('results_merged_absrf_625nm_500x_rep1');
tmp_ab_profile(2,:)=tmp_ab_chlo_results(tmp_ab_chlo_results(:,3)==1,7)/(tmp_rtsum.data(5));
spectrum_retrab(2,1)=tmp_rtsum.data(2)/tmp_rtsum.data(5);%re;
spectrum_retrab(2,2)=tmp_rtsum.data(3)/tmp_rtsum.data(5);%re;
spectrum_retrab(2,3)=tmp_rtsum.data(1)/tmp_rtsum.data(5);%re;

ab_profile=[0.1,0.9]*tmp_ab_profile;
ab_all=sum(ab_profile);

%bar(tmp_ab_profile')
save('ab_profile.mat','ab_profile')
