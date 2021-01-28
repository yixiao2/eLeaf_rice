% eLeaf: 3D model of rice leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 1.2.5

%function e_study_v1_2_3_forfvcb_fit(ARG_TYPE)
tic
%% study simulate measurement A and PhiPSII
if ARG_TYPE==0
    suffix='CKIR64';
elseif ARG_TYPE==1
    suffix='HCIR64';
end
%export_meas4comsol_new_selected(ARG_TYPE);
for loop_std=[1,2]%1:4
    clear tmp_study tmp_std1
    %% study1~4
    %load study1
    %tmp_study=load(['study',num2str(loop_std),'_',suffix]);
    if loop_std==1||loop_std==3
        tmp_study(1,:)=1e-2*ones(1,14);
        tmp_study(2,:)=[25,50,75,150,200,300,400,600,800,1000,1200,1400,1600,2000]*1.0e-6;
    else
        tmp_study(1,:)=[50,80,120,150,230,290,360,500,650,720,890,1150]/10*0.33e-3;
        tmp_study(2,:)=2000.0*1e-6*ones(1,12);
    end
    tmp_std1{1}=num2str(tmp_study(1,:));
    tmp_std1{2}=num2str(tmp_study(2,:));
    if loop_std==1||loop_std==2 %% study1 and study2 are normal O2
        tmp_std1{3}=num2str(ones(size(tmp_study(1,:)))*21000);
    else
        tmp_std1{3}=num2str(ones(size(tmp_study(1,:)))*2100);
    end
    tmpstr_std=['std',num2str(loop_std)];
    model.study.create(tmpstr_std);
    model.study(tmpstr_std).create('stat', 'Stationary');
    
    tmpstr_sol=['sol',num2str(loop_std)];
    model.sol.create(tmpstr_sol);
    model.sol(tmpstr_sol).study(tmpstr_std);
    model.sol(tmpstr_sol).attach(tmpstr_std);
    model.sol(tmpstr_sol).create('st1', 'StudyStep');
    model.sol(tmpstr_sol).create('v1', 'Variables');
    model.sol(tmpstr_sol).create('s1', 'Stationary');
    model.sol(tmpstr_sol).feature('s1').create('fc1', 'FullyCoupled');
    model.sol(tmpstr_sol).feature('s1').create('pDef', 'Parametric');
    
    model.study(tmpstr_std).feature('stat').set('punit', {'mol/m^3' 'mol/m^2/s' 'Pa'});
    model.study(tmpstr_std).feature('stat').set('preusesol', 'yes');
    %model.study('std1').feature('stat').set('plistarr', {'0.0099, 0.0198, 0.0099, 0.0198, 0.0099, 0.0198, 0.0099, 0.0198, 0.0099, 0.0198, 0.0099, 0.0198, 0.0099, 0.0198, 0.0099, 0.0198, 0.0099, 0.0198, 0.0099, 0.0198, 0.0099, 0.0198, 0.0099, 0.0198, 2.0E-4, 0.0020, 0.0040, 0.0060, 0.01, 0.014, 0.016, 0.018, 0.02, 0.024, 0.028, 0.033, 0.04' '50e-6,50e-6,100e-6,100e-6,150e-6,150e-6,200e-6,200e-6,300e-6,300e-6,400e-6,400e-6,600e-6,600e-6,800e-6,800e-6,1000e-6,1000e-6,1300e-6,1300e-6,1600e-6,1600e-6,2000e-6,2000e-6,2000e-6,2000e-6,2000e-6,2000e-6,2000e-6,2000e-6,2000e-6,2000e-6,2000e-6,2000e-6,2000e-6,2000e-6,2000e-6'});
    model.study(tmpstr_std).feature('stat').set('plistarr',tmp_std1);
    model.study(tmpstr_std).feature('stat').set('useparam', 'on');
    model.study(tmpstr_std).feature('stat').set('pcontinuationmode', 'no');
    model.study(tmpstr_std).feature('stat').set('pname', {'Ci' 'irra' 'pO2'});
    
    model.sol(tmpstr_sol).attach(tmpstr_std);
    %model.sol('sol1').feature('v1').set('clist', {'0.0099, 0.0198, 0.0099, 0.0198, 0.0099, 0.0198, 0.0099, 0.0198, 0.0099, 0.0198, 0.0099, 0.0198, 0.0099, 0.0198, 0.0099, 0.0198, 0.0099, 0.0198, 0.0099, 0.0198, 0.0099, 0.0198, 0.0099, 0.0198, 2.0E-4, 0.0020, 0.0040, 0.0060, 0.01, 0.014, 0.016, 0.018, 0.02, 0.024, 0.028, 0.033, 0.04' '50e-6,50e-6,100e-6,100e-6,150e-6,150e-6,200e-6,200e-6,300e-6,300e-6,400e-6,400e-6,600e-6,600e-6,800e-6,800e-6,1000e-6,1000e-6,1300e-6,1300e-6,1600e-6,1600e-6,2000e-6,2000e-6,2000e-6,2000e-6,2000e-6,2000e-6,2000e-6,2000e-6,2000e-6,2000e-6,2000e-6,2000e-6,2000e-6,2000e-6,2000e-6'});
    model.sol(tmpstr_sol).feature('v1').set('clist', tmp_std1);
    model.sol(tmpstr_sol).feature('v1').set('cname', {'Ci' 'irra' 'pO2'});
    model.sol(tmpstr_sol).feature('v1').set('clistctrl', {'pDef' 'pDef'});
    model.sol(tmpstr_sol).feature('s1').set('probesel', 'none');
    model.sol(tmpstr_sol).feature('s1').feature('dDef').set('linsolver', 'pardiso');
    model.sol(tmpstr_sol).feature('s1').feature('fc1').set('maxiter', '40');
    model.sol(tmpstr_sol).feature('s1').feature('pDef').set('uselsqdata', 'off');
    model.sol(tmpstr_sol).feature('s1').feature('pDef').set('ponerror', 'empty');
    model.sol(tmpstr_sol).feature('s1').feature('pDef').set('preusesol', 'yes');
    model.sol(tmpstr_sol).feature('s1').feature('pDef').set('pname', {'Ci' 'irra' 'pO2'});
    model.sol(tmpstr_sol).feature('s1').feature('pDef').set('punit', {'mol/m^3' 'mol/m^2/s' 'Pa'});
    %model.sol('sol1').feature('s1').feature('pDef').set('plistarr', {'0.0099, 0.0198, 0.0099, 0.0198, 0.0099, 0.0198, 0.0099, 0.0198, 0.0099, 0.0198, 0.0099, 0.0198, 0.0099, 0.0198, 0.0099, 0.0198, 0.0099, 0.0198, 0.0099, 0.0198, 0.0099, 0.0198, 0.0099, 0.0198, 2.0E-4, 0.0020, 0.0040, 0.0060, 0.01, 0.014, 0.016, 0.018, 0.02, 0.024, 0.028, 0.033, 0.04' '50e-6,50e-6,100e-6,100e-6,150e-6,150e-6,200e-6,200e-6,300e-6,300e-6,400e-6,400e-6,600e-6,600e-6,800e-6,800e-6,1000e-6,1000e-6,1300e-6,1300e-6,1600e-6,1600e-6,2000e-6,2000e-6,2000e-6,2000e-6,2000e-6,2000e-6,2000e-6,2000e-6,2000e-6,2000e-6,2000e-6,2000e-6,2000e-6,2000e-6,2000e-6'});
    model.sol(tmpstr_sol).feature('s1').feature('pDef').set('plistarr', tmp_std1);
    model.sol(tmpstr_sol).feature('s1').feature('pDef').set('pcontinuationmode', 'no');
    model.sol(tmpstr_sol).runAll;
end

%% output
count_tbl=1;
count_int=1;
count_av=1;
for loop_std=[1,2]%1:4
    % A of std1
    tmp_str_int=['int',num2str(count_int)];
    count_int=count_int+1;
    model.result.numerical.create(tmp_str_int, 'IntVolume');
    
    model.result.numerical(tmp_str_int).selection.named('uni1');
    model.result.numerical(tmp_str_int).setIndex('expr', '(fvc-integ_src_fphres)/leafsurface', 0);
    model.result.numerical(tmp_str_int).set('data', ['dset',num2str(loop_std)]);
    
    tmp_str1=['tbl',num2str(count_tbl)];
    model.result.table.create(tmp_str1, 'Table');
    count_tbl=count_tbl+1;
    model.result.table(tmp_str1).comments('Volume Integration 1 ((fvc-integ_src_fphres)/leafsurface)');
    model.result.numerical(tmp_str_int).set('table', tmp_str1);
    model.result.numerical(tmp_str_int).setResult;
    tmp_str2=['a_std',num2str(loop_std),'_prswp_',suffix,'.txt'];%'a_std1.txt'
    model.result.table(tmp_str1).save(tmp_str2);
    
    % PhiPSII of std1
    tmp_str2='bschlo';
    tmp_str_int=['int',num2str(count_int)];
    count_int=count_int+1;
    model.result.numerical.create(tmp_str_int, 'IntVolume');
    model.result.numerical(tmp_str_int).selection.named(tmp_str2);
    model.result.numerical(tmp_str_int).setIndex('expr', '(fvc)/leafsurface', 0);
    model.result.numerical(tmp_str_int).set('data', ['dset',num2str(loop_std)]);
    
    tmp_str1=['tbl',num2str(count_tbl)];
    model.result.table.create(tmp_str1, 'Table');
    count_tbl=count_tbl+1;
    model.result.table(tmp_str1).comments('fcv/leafsurface and Cc');
    model.result.numerical(tmp_str_int).set('table', tmp_str1);
    model.result.numerical(tmp_str_int).setResult;
    
    tmp_str_av=['av',num2str(count_av)];
    count_av=count_av+1;
    model.result.numerical.create(tmp_str_av, 'AvVolume');
    model.result.numerical(tmp_str_av).set('data', ['dset',num2str(loop_std)]);
    model.result.numerical(tmp_str_av).selection.named(tmp_str2);
    model.result.numerical(tmp_str_av).setIndex('expr', 'c1/sc', 0);
    model.result.numerical(tmp_str_av).set('table', tmp_str1);
    model.result.numerical(tmp_str_av).appendResult;
    
    for tmp_loop=1:count_chl
        tmp_str2=['mschlo_',num2str(tmp_loop)];
        %model.result.numerical.create('int2', 'IntVolume');
        tmp_str_int=['int',num2str(count_int)];
        count_int=count_int+1;
        model.result.numerical.create(tmp_str_int, 'IntVolume');
        model.result.numerical(tmp_str_int).selection.named(tmp_str2);
        model.result.numerical(tmp_str_int).setIndex('expr', '(fvc)/leafsurface', 0);
        model.result.numerical(tmp_str_int).set('data', ['dset',num2str(loop_std)]);
        %model.result.table.create('tbl2', 'Table');
        model.result.table(tmp_str1).comments('fcv/leafsurface and Cc');
        model.result.numerical(tmp_str_int).set('table', tmp_str1);
        model.result.numerical(tmp_str_int).appendResult;
        %model.result.numerical.create('av1', 'AvVolume');
        tmp_str_av=['av',num2str(count_av)];
        count_av=count_av+1;
        model.result.numerical.create(tmp_str_av, 'AvVolume');
        model.result.numerical(tmp_str_av).set('data', ['dset',num2str(loop_std)]);
        model.result.numerical(tmp_str_av).selection.named(tmp_str2);
        model.result.numerical(tmp_str_av).setIndex('expr', 'c1/sc', 0);
        model.result.numerical(tmp_str_av).set('table', tmp_str1);
        model.result.numerical(tmp_str_av).appendResult;
    end
    tmp_str2=['phipsii_std',num2str(loop_std),'_prswp_',suffix,'.txt'];%'phipsii_std1.txt'
    model.result.table(tmp_str1).save(tmp_str2);

    tmp_str2='bschlo';
    tmp_str1=['tbl',num2str(count_tbl)];
    model.result.table.create(tmp_str1, 'Table');
    count_tbl=count_tbl+1;
    %%cc
    tmp_str_av=['av',num2str(count_av)];
    count_av=count_av+1;
    model.result.numerical.create(tmp_str_av, 'AvVolume');
    model.result.numerical(tmp_str_av).set('data', ['dset',num2str(loop_std)]);
    model.result.numerical(tmp_str_av).selection.named(tmp_str2);
    model.result.numerical(tmp_str_av).setIndex('expr', 'c1/sc', 0);
    model.result.numerical(tmp_str_av).set('table', tmp_str1);
    model.result.numerical(tmp_str_av).setResult;
    %%fc
    tmp_str_int=['int',num2str(count_int)];
    count_int=count_int+1;
    model.result.numerical.create(tmp_str_int, 'IntVolume');
    model.result.numerical(tmp_str_int).selection.named(tmp_str2);
    model.result.numerical(tmp_str_int).setIndex('expr', '(fc)/leafsurface', 0);
    model.result.numerical(tmp_str_int).set('data', ['dset',num2str(loop_std)]);
    model.result.numerical(tmp_str_int).set('table', tmp_str1);
    model.result.numerical(tmp_str_int).appendResult;
    %%fj
    tmp_str_int=['int',num2str(count_int)];
    count_int=count_int+1;
    model.result.numerical.create(tmp_str_int, 'IntVolume');
    model.result.numerical(tmp_str_int).selection.named(tmp_str2);
    model.result.numerical(tmp_str_int).setIndex('expr', '(fj)/leafsurface', 0);
    model.result.numerical(tmp_str_int).set('data', ['dset',num2str(loop_std)]);
    model.result.numerical(tmp_str_int).set('table', tmp_str1);
    model.result.numerical(tmp_str_int).appendResult;
    %%vm
    tmp_str_int=['int',num2str(count_int)];
    count_int=count_int+1;
    model.result.numerical.create(tmp_str_int, 'IntVolume');
    model.result.numerical(tmp_str_int).selection.named(tmp_str2);
    model.result.numerical(tmp_str_int).setIndex('expr', '(Vcmax_fvcb/Vols_all)', 0);
    model.result.numerical(tmp_str_int).set('data', ['dset',num2str(loop_std)]);
    model.result.numerical(tmp_str_int).set('table', tmp_str1);
    model.result.numerical(tmp_str_int).appendResult;
    %%j
    tmp_str_int=['int',num2str(count_int)];
    count_int=count_int+1;
    model.result.numerical.create(tmp_str_int, 'IntVolume');
    model.result.numerical(tmp_str_int).selection.named(tmp_str2);
    model.result.numerical(tmp_str_int).setIndex('expr', '(J)/leafsurface', 0);
    model.result.numerical(tmp_str_int).set('data', ['dset',num2str(loop_std)]);
    model.result.numerical(tmp_str_int).set('table', tmp_str1);
    model.result.numerical(tmp_str_int).appendResult;
    
    for tmp_loop=1:count_chl
        tmp_str2=['mschlo_',num2str(tmp_loop)];
        %tmp_str1=['tbl',num2str(count_tbl)];
        %model.result.table.create(tmp_str1, 'Table');
        %count_tbl=count_tbl+1;
        %%cc
        tmp_str_av=['av',num2str(count_av)];
        count_av=count_av+1;
        model.result.numerical.create(tmp_str_av, 'AvVolume');
        model.result.numerical(tmp_str_av).set('data', ['dset',num2str(loop_std)]);
        model.result.numerical(tmp_str_av).selection.named(tmp_str2);
        model.result.numerical(tmp_str_av).setIndex('expr', 'c1/sc', 0);
        model.result.numerical(tmp_str_av).set('table', tmp_str1);
        model.result.numerical(tmp_str_av).appendResult;
        %%fc
        tmp_str_int=['int',num2str(count_int)];
        count_int=count_int+1;
        model.result.numerical.create(tmp_str_int, 'IntVolume');
        model.result.numerical(tmp_str_int).selection.named(tmp_str2);
        model.result.numerical(tmp_str_int).setIndex('expr', '(fc)/leafsurface', 0);
        model.result.numerical(tmp_str_int).set('data', ['dset',num2str(loop_std)]);
        model.result.numerical(tmp_str_int).set('table', tmp_str1);
        model.result.numerical(tmp_str_int).appendResult;
        %%fj
        tmp_str_int=['int',num2str(count_int)];
        count_int=count_int+1;
        model.result.numerical.create(tmp_str_int, 'IntVolume');
        model.result.numerical(tmp_str_int).selection.named(tmp_str2);
        model.result.numerical(tmp_str_int).setIndex('expr', '(fj)/leafsurface', 0);
        model.result.numerical(tmp_str_int).set('data', ['dset',num2str(loop_std)]);
        model.result.numerical(tmp_str_int).set('table', tmp_str1);
        model.result.numerical(tmp_str_int).appendResult;
        %%vm
        tmp_str_int=['int',num2str(count_int)];
        count_int=count_int+1;
        model.result.numerical.create(tmp_str_int, 'IntVolume');
        model.result.numerical(tmp_str_int).selection.named(tmp_str2);
        model.result.numerical(tmp_str_int).setIndex('expr', '(Vcmax_fvcb/Vols_all)', 0);
        model.result.numerical(tmp_str_int).set('data', ['dset',num2str(loop_std)]);
        model.result.numerical(tmp_str_int).set('table', tmp_str1);
        model.result.numerical(tmp_str_int).appendResult;
        %%j
        tmp_str_int=['int',num2str(count_int)];
        count_int=count_int+1;
        model.result.numerical.create(tmp_str_int, 'IntVolume');
        model.result.numerical(tmp_str_int).selection.named(tmp_str2);
        model.result.numerical(tmp_str_int).setIndex('expr', '(J)/leafsurface', 0);
        model.result.numerical(tmp_str_int).set('data', ['dset',num2str(loop_std)]);
        model.result.numerical(tmp_str_int).set('table', tmp_str1);
        model.result.numerical(tmp_str_int).appendResult;
    end
    tmp_str2=['cov_std',num2str(loop_std),'_prswp_',suffix,'.txt'];%'cov_std1.txt'
    model.result.table(tmp_str1).save(tmp_str2);
end

mphsave(model,['eleaf_fvcb_prswp_',suffix,'_a_phipsii.mph'])
toc
%end
