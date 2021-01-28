% eLeaf: 3D model of rice leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 1.2.5

%% process dissection; export txt for extended FvCB model
function recalcov_v_1_2_3(tmp_suffix)
model=mphload(['eleaf_fvcb_prswp_',tmp_suffix,'_a_phipsii.mph']);
load save_e_geom.mat count_chl
%%get count_int count_av
tmp_tags=model.result.numerical.tags;
count_int=0;
count_av=0;
loop=1;
if size(tmp_tags,1)~=0
    while count_int==0||count_av==0||loop~=size(tmp_tags,1)
        tmp_str=char(tmp_tags(end-loop+1));
        if isempty(regexp(tmp_str,'^int\d+$','match'))==0&&count_int==0
            tmp_str2=regexp(tmp_str,'\d+$','match');
            tmp_str2=tmp_str2{1};
            count_int=str2num(tmp_str2);
        elseif isempty(regexp(tmp_str,'^av\d+$','match'))==0&&count_av==0
            tmp_str2=regexp(tmp_str,'\d+$','match');
            tmp_str2=tmp_str2{1};
            count_av=str2num(tmp_str2);
        end
        loop=loop+1;
    end
end
count_int=count_int+1;
count_av=count_av+1;
%%get count_tbl
count_tbl=0;
tmp_tags=model.result.table.tags;
loop=1;
if size(tmp_tags,1)~=1
    while count_tbl==0||loop~=size(tmp_tags,1)
        tmp_str=char(tmp_tags(end-loop+1));
        if isempty(regexp(tmp_str,'^tbl\d+$','match'))==0&&count_tbl==0
            tmp_str2=regexp(tmp_str,'\d+$','match');
            tmp_str2=tmp_str2{1};
            count_tbl=str2num(tmp_str2);
        end
        loop=loop+1;
    end
end
count_tbl=count_tbl+1;

%% new calculation and new tables for dissect coordination
%% need to export {cc,fc,fj,vm,j}
for loop_std=1:4
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
    tmp_str2=['cov_std',num2str(loop_std),'_prswp_',tmp_suffix,'.txt'];%'cov_std1.txt'
    model.result.table(tmp_str1).save(tmp_str2);
end
end
