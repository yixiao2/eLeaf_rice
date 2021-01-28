% eLeaf: 3D model of rice leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 1.2.4

function [tmp_wallString,tmp_chloString,tmp_chliString,tmp_IDXxse,tmp_IDXyse]=msoutline_v1_2(tmp_Nlobe,dw,dse,epsln)
tmp_Tlobe=tmp_Nlobe/2;
tmp_Npts=14+4*(tmp_Tlobe-2)+2*tmp_Nlobe;
%#coordinate x of points
tmp_count=1;%point 1
tmp_IDXx{tmp_count}=['1/2*(',num2str(-tmp_Tlobe),')'];
tmp_IDXy{tmp_count}=['1/(2*sqrt(3))*(',num2str(1),')'];%1/(2*sqrt(3))*(1)
tmp_IDXxw{tmp_count}=['(1/2*(-',num2str(tmp_Tlobe),'))*(1-2*dw)+(1/2*(-',num2str(tmp_Tlobe),'+1))*(2*dw)'];%[1]1/2*(-Tlobe);[2]1/2*(-Tlobe+1) under 2*dw
tmp_IDXyw{tmp_count}=['(1/(2*sqrt(3)))*(1-dw*2)'];%[1]1/sqrt(3);[2]0 under 2*dw
if tmp_Nlobe==4
    if dse<(1/2-dw)/(1-dw+epsln/2)
        tmp_IDXxse{tmp_count}=['(1/2*(-',num2str(tmp_Tlobe),')+dw)*(1-dse)+(1/2*(-',num2str(tmp_Tlobe),'+2)+epsln4/2)*(dse)'];%[1]1/2*(-Tlobe)+dw;[2]1/2*(-Tlobe+2)+epsln4/2 under dse
        tmp_IDXyse{tmp_count}=['(1/(2*sqrt(3)))*(1-(dw+dse*(1-dw+epsln4/2))/(1/2))'];%[1]1/(2*sqrt(3));[2]0 under (dw+dse*(1-dw+epsln4/2))/(1/2)
    else
        tmp_IDXxse{tmp_count}=['(1/2*(-',num2str(tmp_Tlobe),')+dw)*(1-dse)+(1/2*(-',num2str(tmp_Tlobe),'+2)+epsln4/2)*(dse)'];%[1]1/2*(-Tlobe)+dw;[2]1/2*(-Tlobe+2)+epsln4/2 under dse
        tmp_IDXyse{tmp_count}=['0'];
    end
else
    if dse<(1/2-dw)/(1-dw+epsln/2)
        tmp_IDXxse{tmp_count}=['(1/2*(-',num2str(tmp_Tlobe),')+dw)*(1-dse)+(1/2*(-',num2str(tmp_Tlobe),'+2)+epsln/2)*(dse)'];%[1]1/2*(-Tlobe)+dw;[2]1/2*(-Tlobe+2)+epsln/2 under dse
        tmp_IDXyse{tmp_count}=['(1/(2*sqrt(3)))*(1-(dw+dse*(1-dw+epsln/2))/(1/2))'];%[1]1/(2*sqrt(3));[2]0 under (dw+dse*(1-dw+epsln/2))/(1/2)
    else
        tmp_IDXxse{tmp_count}=['(1/2*(-',num2str(tmp_Tlobe),')+dw)*(1-dse)+(1/2*(-',num2str(tmp_Tlobe),'+2)+epsln/2)*(dse)'];%[1]1/2*(-Tlobe)+dw;[2]1/2*(-Tlobe+2)+epsln/2 under dse
        tmp_IDXyse{tmp_count}=['0'];
    end
end
tmp_count=tmp_count+1;%end point 1
%point 2
tmp_IDXx{tmp_count}=['1/2*(',num2str(-tmp_Tlobe),')*(lmbd)+1/2*(',num2str(-tmp_Tlobe+1),')*(1-lmbd)'];%1/2*(-4)*(lmbd)+1/2*(-3)*(1-lmbd)
tmp_IDXy{tmp_count}=['1/(2*sqrt(3))*(',num2str(1),')*(lmbd)+1/(2*sqrt(3))*(',num2str(2),')*(1-lmbd)'];%1/(2*sqrt(3))*(1)*(lmbd)+1/(2*sqrt(3))*(2)*(1-lmbd)
tmp_IDXxw{tmp_count}=['(1/2*(',num2str(-tmp_Tlobe),')*(lmbd)+1/2*(',num2str(-tmp_Tlobe+1),')*(1-lmbd))*(1-2*dw)+(1/2*(',num2str(-tmp_Tlobe+1),'))*(2*dw)'];%[1]1/2*(-Tlobe)*(lmbd)+1/2*(-Tlobe+1)*(1-lmbd);[2]1/2*(-Tlobe+1) under 2*dw
tmp_IDXyw{tmp_count}=['(1/(2*sqrt(3))*(1)*(lmbd)+1/(2*sqrt(3))*(2)*(1-lmbd))*(1-2*dw)'];%[1]1/(2*sqrt(3))*(1)*(lmbd)+1/(2*sqrt(3))*(2)*(1-lmbd);[2]1/(2*sqrt(3)) under 2*dw
if tmp_Nlobe==4
    if dse<(1/2-dw)/(1-dw+epsln/2)
        tmp_IDXxse{tmp_count}=['(1/2*(-',num2str(tmp_Tlobe),')*(lmbd)+1/2*(-',num2str(tmp_Tlobe),'+1)*(1-lmbd))*(1-dse*(1-dw+epsln4/2)/(1/2-dw))+(1/2*(-',num2str(tmp_Tlobe),'+1))*(dse*(1-dw+epsln4/2)/(1/2-dw))'];%[1]1/2*(-Tlobe)*(lmbd)+1/2*(-Tlobe+1)*(1-lmbd);[2]1/2*(-Tlobe+1) under dse*(1-dw+epsln4/2)/(1/2-dw)
        tmp_IDXyse{tmp_count}=['(1/(2*sqrt(3))*(1)*(lmbd)+1/(2*sqrt(3))*(2)*(1-lmbd))*(1-(dw+dse*(1-dw+epsln4/2))/(1/2))'];%[1]1/(2*sqrt(3))*(1)*(lmbd)+1/(2*sqrt(3))*(2)*(1-lmbd);[2]0 under (dw+dse*(1-dw+epsln4/2))/(1/2)
    else
        tmp_IDXxse{tmp_count}=tmp_IDXxse{tmp_count-1};
        tmp_IDXyse{tmp_count}=['0'];
    end
else
    if dse<(1/2-dw)/(1-dw+epsln/2)
        tmp_IDXxse{tmp_count}=['(1/2*(-',num2str(tmp_Tlobe),')*(lmbd)+1/2*(-',num2str(tmp_Tlobe),'+1)*(1-lmbd))*(1-dse*(1-dw+epsln/2)/(1/2-dw))+(1/2*(-',num2str(tmp_Tlobe),'+1))*(dse*(1-dw+epsln/2)/(1/2-dw))'];%[1]1/2*(-Tlobe)*(lmbd)+1/2*(-Tlobe+1)*(1-lmbd);[2]1/2*(-Tlobe+1) under dse*(1-dw+epsln/2)/(1/2-dw)
        tmp_IDXyse{tmp_count}=['(1/(2*sqrt(3))*(1)*(lmbd)+1/(2*sqrt(3))*(2)*(1-lmbd))*(1-(dw+dse*(1-dw+epsln/2))/(1/2))'];%[1]1/(2*sqrt(3))*(1)*(lmbd)+1/(2*sqrt(3))*(2)*(1-lmbd);[2]0 under (dw+dse*(1-dw+epsln/2))/(1/2)
    else
        tmp_IDXxse{tmp_count}=tmp_IDXxse{tmp_count-1};
        tmp_IDXyse{tmp_count}=['0'];
    end
end
tmp_count=tmp_count+1;%end point 2
%point 3
tmp_IDXx{tmp_count}=['1/2*(',num2str(-tmp_Tlobe+2),')*(rho)+1/2*(',num2str(-tmp_Tlobe+1),')*(1-rho)'];%1/2*(-2)*(rho)+1/2*(-3)*(1-rho)
tmp_IDXy{tmp_count}=['1/(2*sqrt(3))*(',num2str(1),')*(rho)+1/(2*sqrt(3))*(',num2str(2),')*(1-rho)'];%1/(2*sqrt(3))*(1)*(rho)+1/(2*sqrt(3))*(2)*(1-rho)
tmp_IDXxw{tmp_count}=['(1/2*(-',num2str(tmp_Tlobe),'+2)*(rho)+1/2*(-',num2str(tmp_Tlobe),'+1)*(1-rho))*(1-dw/(1/2*(1-rho)))+(1/2*(-',num2str(tmp_Tlobe),'+2))*(dw/(1/2*(1-rho)))'];%[1]1/2*(-Tlobe+2)*(rho)+1/2*(-Tlobe+1)*(1-rho);[2]1/2*(-Tlobe+2) under dw/(1/2*(1-rho))
tmp_IDXyw{tmp_count}=['(1/(2*sqrt(3))*(1)*(rho)+1/(2*sqrt(3))*(2)*(1-rho))*(1-dw/(1/2*(1-rho)))+(1/(2*sqrt(3)))*(dw/(1/2*(1-rho)))'];%[1]1/(2*sqrt(3))*(1)*(rho)+1/(2*sqrt(3))*(2)*(1-rho);[2]1/(2*sqrt(3)) under dw/(1/2*(1-rho))
if tmp_Nlobe==4
    tmp_IDXxse{tmp_count}=['(1/2*(-',num2str(tmp_Tlobe),'+2)*(rho)+1/2*(-',num2str(tmp_Tlobe),'+1)*(1-rho))*(1-((1/sqrt(3)*(1-rho)+epsln4/sqrt(3)-2*dw)*dse+2*dw)/(1/sqrt(3)*(1-rho)))+(1/2*(-',num2str(tmp_Tlobe),'+2))*(((1/sqrt(3)*(1-rho)+epsln4/sqrt(3)-2*dw)*dse+2*dw)/(1/sqrt(3)*(1-rho)))'];%[1]1/2*(-Tlobe+2)*(rho)+1/2*(-Tlobe+1)*(1-rho);[2]1/2*(-Tlobe+2) under ((1/sqrt(3)*(1-rho)+epsln4/sqrt(3)-2*dw)*dse+2*dw)/(1/sqrt(3)*(1-rho))
    tmp_IDXyse{tmp_count}=['(1/(2*sqrt(3))*(1)*(rho)+1/(2*sqrt(3))*(2)*(1-rho))*(1-((1/sqrt(3)*(1-rho)+epsln4/sqrt(3)-2*dw)*dse+2*dw)/(1/sqrt(3)*(1-rho)))+(1/(2*sqrt(3)))*(((1/sqrt(3)*(1-rho)+epsln4/sqrt(3)-2*dw)*dse+2*dw)/(1/sqrt(3)*(1-rho)))'];%[1]1/(2*sqrt(3))*(1)*(rho)+1/(2*sqrt(3))*(2)*(1-rho);[2]1/(2*sqrt(3)) under ((1/sqrt(3)*(1-rho)+epsln4/sqrt(3)-2*dw)*dse+2*dw)/(1/sqrt(3)*(1-rho))
else
    tmp_IDXxse{tmp_count}=['(1/2*(-',num2str(tmp_Tlobe),'+2)*(rho)+1/2*(-',num2str(tmp_Tlobe),'+1)*(1-rho))*(1-((1/sqrt(3)*(1-rho)+epsln/sqrt(3)-2*dw)*dse+2*dw)/(1/sqrt(3)*(1-rho)))+(1/2*(-',num2str(tmp_Tlobe),'+2))*(((1/sqrt(3)*(1-rho)+epsln/sqrt(3)-2*dw)*dse+2*dw)/(1/sqrt(3)*(1-rho)))'];%[1]1/2*(-Tlobe+2)*(rho)+1/2*(-Tlobe+1)*(1-rho);[2]1/2*(-Tlobe+2) under ((1/sqrt(3)*(1-rho)+epsln/sqrt(3)-2*dw)*dse+2*dw)/(1/sqrt(3)*(1-rho))
    tmp_IDXyse{tmp_count}=['(1/(2*sqrt(3))*(1)*(rho)+1/(2*sqrt(3))*(2)*(1-rho))*(1-((1/sqrt(3)*(1-rho)+epsln/sqrt(3)-2*dw)*dse+2*dw)/(1/sqrt(3)*(1-rho)))+(1/(2*sqrt(3)))*(((1/sqrt(3)*(1-rho)+epsln/sqrt(3)-2*dw)*dse+2*dw)/(1/sqrt(3)*(1-rho)))'];%[1]1/(2*sqrt(3))*(1)*(rho)+1/(2*sqrt(3))*(2)*(1-rho);[2]1/(2*sqrt(3)) under ((1/sqrt(3)*(1-rho)+epsln/sqrt(3)-2*dw)*dse+2*dw)/(1/sqrt(3)*(1-rho))
end
tmp_count=tmp_count+1;%end point3
if tmp_Nlobe==4
    %point 4
    tmp_IDXx{tmp_count}=['1/2*(',num2str(-tmp_Tlobe+1),')'];
    tmp_IDXy{tmp_count}=['1/(2*sqrt(3))*(',num2str(4),')*(lmbd)+1/(2*sqrt(3))*(',num2str(2),')*(1-lmbd)'];%1/(2*sqrt(3))*(4)*(lmbd)+1/(2*sqrt(3))*(2)*(1-lmbd)
    tmp_IDXxw{tmp_count}=['(1/2*(',num2str(-tmp_Tlobe+1),'))*(1-2*dw)+(1/2*(',num2str(-tmp_Tlobe+2),'))*(2*dw)'];%[1]1/2*(-Tlobe+1);[2]1/2*(-Tlobe+2) under 2*dw
    tmp_IDXyw{tmp_count}=['(1/(2*sqrt(3))*(4)*(lmbd)+1/(2*sqrt(3))*(2)*(1-lmbd))*(1-2*dw)+(1/(2*sqrt(3))*3)*(2*dw)'];%[1]1/(2*sqrt(3))*(4)*(lmbd)+1/(2*sqrt(3))*(2)*(1-lmbd);[2]1/(2*sqrt(3))*3 under 2*dw
    if dse<(1/2-dw)/(1-dw+epsln/2)
        tmp_IDXxse{tmp_count}=['(1/2*(-',num2str(tmp_Tlobe),'+1))*(1-dse*(1-dw+epsln4/2)/(1/2-dw))+(1/2*(-',num2str(tmp_Tlobe),'+2))*(dse*(1-dw+epsln4/2)/(1/2-dw))'];%[1]1/2*(-Tlobe+1);[2]1/2*(-Tlobe+2) under dse*(1-dw+epsln4/2)/(1/2-dw)
        tmp_IDXyse{tmp_count}=['(1/(2*sqrt(3))*(4)*(lmbd)+1/(2*sqrt(3))*(2)*(1-lmbd))*(1-(dw+dse*(1-dw+epsln4/2))/(1/2))+(1/(2*sqrt(3))*3)*((dw+dse*(1-dw+epsln4/2))/(1/2))'];%[1]1/(2*sqrt(3))*(4)*(lmbd)+1/(2*sqrt(3))*(2)*(1-lmbd);[2]1/(2*sqrt(3))*3 under (dw+dse*(1-dw+epsln4/2))/(1/2)
    else
        tmp_IDXxse{tmp_count}=['1/2*(-',num2str(tmp_Tlobe),'+2)'];
        tmp_IDXyse{tmp_count}=['(1/(2*sqrt(3))*3)*(1-(dse-((1/2-dw)/(1-dw+epsln4/2)))/(1-((1/2-dw)/(1-dw+epsln4/2))))+(1/(2*sqrt(3))*1)*((dse-((1/2-dw)/(1-dw+epsln4/2)))/(1-((1/2-dw)/(1-dw+epsln4/2))))'];%when dse=(1/2-dw)/(1-dw+epsln4/2) yse=1/(2*sqrt(3))*3; when dse=1 yse=1/(2*sqrt(3))*1
    end
    tmp_count=tmp_count+1;%end point 4
    %point 5
    tmp_IDXx{tmp_count}=['1/2*(',num2str(-tmp_Tlobe+1),')'];
    tmp_IDXy{tmp_count}=['1/(2*sqrt(3))*(',num2str(4),')'];
    tmp_IDXxw{tmp_count}=['(1/2*(',num2str(-tmp_Tlobe+1),'))*(1-2*dw)+(1/2*(',num2str(-tmp_Tlobe+2),'))*(2*dw)'];%[1]1/2*(-Tlobe+1);[2]1/2*(-Tlobe+2) under 2*dw
    tmp_IDXyw{tmp_count}=['(1/(2*sqrt(3))*4)*(1-2*dw)+(1/(2*sqrt(3))*3)*(2*dw)'];%[1]1/(2*sqrt(3))*(4)*(lmbd)+1/(2*sqrt(3))*(2)*(1-lmbd);[2]1/(2*sqrt(3))*3 under 2*dw
    if dse<(1/2-dw)/(1-dw+epsln/2)
        tmp_IDXxse{tmp_count}=['(1/2*(-',num2str(tmp_Tlobe),'+1))*(1-dse*(1-dw+epsln4/2)/(1/2-dw))+(1/2*(-',num2str(tmp_Tlobe),'+2))*(dse*(1-dw+epsln4/2)/(1/2-dw))'];%[1]1/2*(-Tlobe+1);[2]1/2*(-Tlobe+2) under dse*(1-dw+epsln4/2)/(1/2-dw)
        tmp_IDXyse{tmp_count}=['(1/(2*sqrt(3))*4)*(1-(dw+dse*(1-dw+epsln4/2))/(1/2))+(1/(2*sqrt(3))*3)*((dw+dse*(1-dw+epsln4/2))/(1/2))'];%[1]1/(2*sqrt(3))*4;[2]1/(2*sqrt(3))*3 under (dw+dse*(1-dw+epsln4/2))/(1/2)
    else
        tmp_IDXxse{tmp_count}=['1/2*(-',num2str(tmp_Tlobe),'+2)'];
        tmp_IDXyse{tmp_count}=['(1/(2*sqrt(3))*3)*(1-(dse-((1/2-dw)/(1-dw+epsln4/2)))/(1-((1/2-dw)/(1-dw+epsln4/2))))+(1/(2*sqrt(3))*1)*((dse-((1/2-dw)/(1-dw+epsln4/2)))/(1-((1/2-dw)/(1-dw+epsln4/2))))'];%when dse=(1/2-dw)/(1-dw+epsln4/2) yse=1/(2*sqrt(3))*3; when dse=1 yse=1/(2*sqrt(3))*1
    end
    tmp_count=tmp_count+1;%end point 5
    %point 6
    tmp_IDXx{tmp_count}=['1/2*(',num2str(-tmp_Tlobe+2),')'];
    tmp_IDXy{tmp_count}=['1/(2*sqrt(3))*(',num2str(5),')'];
    tmp_IDXxw{tmp_count}=['1/2*(',num2str(-tmp_Tlobe+2),')'];
    tmp_IDXyw{tmp_count}=['(1/(2*sqrt(3))*5)*(1-2*dw)+(1/(2*sqrt(3))*3)*(2*dw)'];%[1]1/(2*sqrt(3))*(5);[2]1/(2*sqrt(3))*4 under 2*dw
    if dse<(1/2-dw)/(1-dw+epsln/2)
        tmp_IDXxse{tmp_count}=['1/2*(',num2str(-tmp_Tlobe+2),')'];
        tmp_IDXyse{tmp_count}=['(1/(2*sqrt(3))*5)*(1-(dw+dse*(1-dw+epsln4/2))/(1/2))+(1/(2*sqrt(3))*3)*((dw+dse*(1-dw+epsln4/2))/(1/2))'];%[1]1/(2*sqrt(3))*5;[2]1/(2*sqrt(3))*3 under (dw+dse*(1-dw+epsln4/2))/(1/2)
    else
        tmp_IDXxse{tmp_count}=['1/2*(',num2str(-tmp_Tlobe+2),')'];
        tmp_IDXyse{tmp_count}=['(1/(2*sqrt(3))*3)*(1-(dse-((1/2-dw)/(1-dw+epsln4/2)))/(1-((1/2-dw)/(1-dw+epsln4/2))))+(1/(2*sqrt(3))*1)*((dse-((1/2-dw)/(1-dw+epsln4/2)))/(1-((1/2-dw)/(1-dw+epsln4/2))))'];%when dse=(1/2-dw)/(1-dw+epsln4/2) yse=1/(2*sqrt(3))*3; when dse=1 yse=1/(2*sqrt(3))*1
    end
    tmp_count=tmp_count+1;%end point 6
else
    %point 4
    tmp_IDXx{tmp_count}=['1/2*(',num2str(-tmp_Tlobe+1),')'];
    tmp_IDXy{tmp_count}=['1/(2*sqrt(3))*(',num2str(4),')*(lmbd)+1/(2*sqrt(3))*(',num2str(2),')*(1-lmbd)'];%1/(2*sqrt(3))*(4)*(lmbd)+1/(2*sqrt(3))*(2)*(1-lmbd)
    tmp_IDXxw{tmp_count}=['(',tmp_IDXxw{2},')-2/4*1*((',tmp_IDXxw{2},')+sqrt(3)*(',tmp_IDXyw{2},')+(-1/2*(-',num2str(tmp_Tlobe),'+3)))'];
    tmp_IDXyw{tmp_count}=['(',tmp_IDXyw{2},')-2/4*sqrt(3)*((',tmp_IDXxw{2},')+sqrt(3)*(',tmp_IDXyw{2},')+(-1/2*(-',num2str(tmp_Tlobe),'+3)))'];
    tmp_IDXxse{tmp_count}=['(',tmp_IDXxse{2},')-2/4*1*((',tmp_IDXxse{2},')+sqrt(3)*(',tmp_IDXyse{2},')+(-1/2*(-',num2str(tmp_Tlobe),'+3)))'];
    tmp_IDXyse{tmp_count}=['(',tmp_IDXyse{2},')-2/4*sqrt(3)*((',tmp_IDXxse{2},')+sqrt(3)*(',tmp_IDXyse{2},')+(-1/2*(-',num2str(tmp_Tlobe),'+3)))'];
    tmp_count=tmp_count+1;%end point 4
    %point 5
    tmp_IDXx{tmp_count}=['1/2*(',num2str(-tmp_Tlobe+1),')'];
    tmp_IDXy{tmp_count}=['1/(2*sqrt(3))*(',num2str(4),')'];%1/(2*sqrt(3))*(4)
    tmp_IDXxw{tmp_count}=['(',tmp_IDXxw{1},')-2/4*1*((',tmp_IDXxw{1},')+sqrt(3)*(',tmp_IDXyw{1},')+(-1/2*(-',num2str(tmp_Tlobe),'+3)))'];
    tmp_IDXyw{tmp_count}=['(',tmp_IDXyw{1},')-2/4*sqrt(3)*((',tmp_IDXxw{1},')+sqrt(3)*(',tmp_IDXyw{1},')+(-1/2*(-',num2str(tmp_Tlobe),'+3)))'];
    tmp_IDXxse{tmp_count}=['(',tmp_IDXxse{1},')-2/4*1*((',tmp_IDXxse{1},')+sqrt(3)*(',tmp_IDXyse{1},')+(-1/2*(-',num2str(tmp_Tlobe),'+3)))'];
    tmp_IDXyse{tmp_count}=['(',tmp_IDXyse{1},')-2/4*sqrt(3)*((',tmp_IDXxse{1},')+sqrt(3)*(',tmp_IDXyse{1},')+(-1/2*(-',num2str(tmp_Tlobe),'+3)))'];
    tmp_count=tmp_count+1;%end point 5
    %point 6
    tmp_IDXx{tmp_count}=['1/2*(',num2str(-tmp_Tlobe+2),')'];
    tmp_IDXy{tmp_count}=['1/(2*sqrt(3))*(',num2str(5),')'];%1/(2*sqrt(3))*(5)
    tmp_IDXxw{tmp_count}=['(',tmp_IDXxw{5},')-2/4*sqrt(3)*(sqrt(3)*(',tmp_IDXxw{5},')+(',tmp_IDXyw{5},')+(-1/2*sqrt(3)*(-',num2str(tmp_Tlobe),'+3)))'];
    tmp_IDXyw{tmp_count}=['(',tmp_IDXyw{5},')-2/4*1*(sqrt(3)*(',tmp_IDXxw{5},')+(',tmp_IDXyw{5},')+(-1/2*sqrt(3)*(-',num2str(tmp_Tlobe),'+3)))'];
    tmp_IDXxse{tmp_count}=['(',tmp_IDXxse{5},')-2/4*sqrt(3)*(sqrt(3)*(',tmp_IDXxse{5},')+(',tmp_IDXyse{5},')+(-1/2*sqrt(3)*(-',num2str(tmp_Tlobe),'+3)))'];
    tmp_IDXyse{tmp_count}=['(',tmp_IDXyse{5},')-2/4*1*(sqrt(3)*(',tmp_IDXxse{5},')+(',tmp_IDXyse{5},')+(-1/2*sqrt(3)*(-',num2str(tmp_Tlobe),'+3)))'];
    tmp_count=tmp_count+1;%end point 6
    %point 7
    tmp_IDXx{tmp_count}=['(',tmp_IDXx{4},')-2/4*sqrt(3)*(sqrt(3)*(',tmp_IDXx{4},')+(',tmp_IDXy{4},')+(-1/2*sqrt(3)*(-',num2str(tmp_Tlobe),'+3)))'];
    tmp_IDXy{tmp_count}=['(',tmp_IDXy{4},')-2/4*1*(sqrt(3)*(',tmp_IDXx{4},')+(',tmp_IDXy{4},')+(-1/2*sqrt(3)*(-',num2str(tmp_Tlobe),'+3)))'];
    tmp_IDXxw{tmp_count}=['(',tmp_IDXxw{4},')-2/4*sqrt(3)*(sqrt(3)*(',tmp_IDXxw{4},')+(',tmp_IDXyw{4},')+(-1/2*sqrt(3)*(-',num2str(tmp_Tlobe),'+3)))'];
    tmp_IDXyw{tmp_count}=['(',tmp_IDXyw{4},')-2/4*1*(sqrt(3)*(',tmp_IDXxw{4},')+(',tmp_IDXyw{4},')+(-1/2*sqrt(3)*(-',num2str(tmp_Tlobe),'+3)))'];
    tmp_IDXxse{tmp_count}=['(',tmp_IDXxse{4},')-2/4*sqrt(3)*(sqrt(3)*(',tmp_IDXxse{4},')+(',tmp_IDXyse{4},')+(-1/2*sqrt(3)*(-',num2str(tmp_Tlobe),'+3)))'];
    tmp_IDXyse{tmp_count}=['(',tmp_IDXyse{4},')-2/4*1*(sqrt(3)*(',tmp_IDXxse{4},')+(',tmp_IDXyse{4},')+(-1/2*sqrt(3)*(-',num2str(tmp_Tlobe),'+3)))'];
    tmp_count=tmp_count+1;%end point 7
    %point 8
    tmp_IDXx{tmp_count}=['(',tmp_IDXx{3},')-2/4*sqrt(3)*(sqrt(3)*(',tmp_IDXx{3},')+(',tmp_IDXy{3},')+(-1/2*sqrt(3)*(-',num2str(tmp_Tlobe),'+3)))'];
    tmp_IDXy{tmp_count}=['(',tmp_IDXy{3},')-2/4*1*(sqrt(3)*(',tmp_IDXx{3},')+(',tmp_IDXy{3},')+(-1/2*sqrt(3)*(-',num2str(tmp_Tlobe),'+3)))'];
    tmp_IDXxw{tmp_count}=['(',tmp_IDXxw{3},')-2/4*sqrt(3)*(sqrt(3)*(',tmp_IDXxw{3},')+(',tmp_IDXyw{3},')+(-1/2*sqrt(3)*(-',num2str(tmp_Tlobe),'+3)))'];
    tmp_IDXyw{tmp_count}=['(',tmp_IDXyw{3},')-2/4*1*(sqrt(3)*(',tmp_IDXxw{3},')+(',tmp_IDXyw{3},')+(-1/2*sqrt(3)*(-',num2str(tmp_Tlobe),'+3)))'];
    tmp_IDXxse{tmp_count}=['(',tmp_IDXxse{3},')-2/4*sqrt(3)*(sqrt(3)*(',tmp_IDXxse{3},')+(',tmp_IDXyse{3},')+(-1/2*sqrt(3)*(-',num2str(tmp_Tlobe),'+3)))'];
    tmp_IDXyse{tmp_count}=['(',tmp_IDXyse{3},')-2/4*1*(sqrt(3)*(',tmp_IDXxse{3},')+(',tmp_IDXyse{3},')+(-1/2*sqrt(3)*(-',num2str(tmp_Tlobe),'+3)))'];
    tmp_count=tmp_count+1;%end point 8
    %point 9~point 2*Tlobe+2
    sym_x=3;
    while tmp_count~=(2*tmp_Tlobe+2+1)
        if mod(sym_x,2)==1
            tmp_IDXx{tmp_count}=['1/2*(-',num2str(tmp_Tlobe),'+',num2str(sym_x),')*2-(',tmp_IDXx{tmp_count-2},')'];
            tmp_IDXy{tmp_count}=tmp_IDXy{tmp_count-2};
            tmp_IDXxw{tmp_count}=['1/2*(-',num2str(tmp_Tlobe),'+',num2str(sym_x),')*2-(',tmp_IDXxw{tmp_count-2},')'];
            tmp_IDXyw{tmp_count}=tmp_IDXyw{tmp_count-2};
            tmp_IDXxse{tmp_count}=['1/2*(-',num2str(tmp_Tlobe),'+',num2str(sym_x),')*2-(',tmp_IDXxse{tmp_count-2},')'];
            tmp_IDXyse{tmp_count}=tmp_IDXyse{tmp_count-2};
            tmp_count=tmp_count+1;
            tmp_IDXx{tmp_count}=['1/2*(-',num2str(tmp_Tlobe),'+',num2str(sym_x+1),')'];
            tmp_IDXy{tmp_count}=['1/(2*sqrt(3))*5'];
            tmp_IDXxw{tmp_count}=tmp_IDXx{tmp_count};
            tmp_IDXyw{tmp_count}=['1/(2*sqrt(3))*5-dw/sqrt(3)*2'];%[1]1/(2*sqrt(3))*5-dw/sqrt(3)*2
            tmp_IDXxse{tmp_count}=tmp_IDXx{tmp_count};
            %         tmp_IDXyse{tmp_count}=['1/(2*sqrt(3))*5-dw/sqrt(3)*2','-dse*(1/(2*sqrt(3))*3+epsln/sqrt(3)-dw/sqrt(3)*2)'];
            if dse<(1/2-dw)/(1-dw+epsln/2)
                tmp_IDXyse{tmp_count}=tmp_IDXyse{6};
            else
                tmp_IDXyse{tmp_count}=['(1/(2*sqrt(3))*3)*(1-(dse-(1/2-dw)/(1-dw+epsln/2))/(1-(1/2-dw)/(1-dw+epsln/2)))+(1/(2*sqrt(3))*2-epsln/sqrt(3))*((dse-(1/2-dw)/(1-dw+epsln/2))/(1-(1/2-dw)/(1-dw+epsln/2)))'];%[1]1/(2*sqrt(3))*3 when dse=(1/2-dw)/(1-dw+epsln/2);[2] 1/(2*sqrt(3))*2-epsln/sqrt(3) when dse=1; so under (dse-(1/2-dw)/(1-dw+epsln/2))/(1-(1/2-dw)/(1-dw+epsln/2))
            end
            tmp_count=tmp_count+1;
        else
            tmp_IDXx{tmp_count}=['1/2*(-',num2str(tmp_Tlobe),'+',num2str(sym_x),')*2-(',tmp_IDXx{tmp_count-2},')'];
            tmp_IDXy{tmp_count}=tmp_IDXy{tmp_count-2};
            tmp_IDXxw{tmp_count}=['1/2*(-',num2str(tmp_Tlobe),'+',num2str(sym_x),')*2-(',tmp_IDXxw{tmp_count-2},')'];
            tmp_IDXyw{tmp_count}=tmp_IDXyw{tmp_count-2};
            tmp_IDXxse{tmp_count}=['1/2*(-',num2str(tmp_Tlobe),'+',num2str(sym_x),')*2-(',tmp_IDXxse{tmp_count-2},')'];
            tmp_IDXyse{tmp_count}=tmp_IDXyse{tmp_count-2};
            tmp_count=tmp_count+1;
            tmp_IDXx{tmp_count}=['1/2*(-',num2str(tmp_Tlobe),'+',num2str(sym_x),')*2-(',tmp_IDXx{tmp_count-4},')'];
            tmp_IDXy{tmp_count}=tmp_IDXy{tmp_count-4};
            tmp_IDXxw{tmp_count}=['1/2*(-',num2str(tmp_Tlobe),'+',num2str(sym_x),')*2-(',tmp_IDXxw{tmp_count-4},')'];
            tmp_IDXyw{tmp_count}=tmp_IDXyw{tmp_count-4};
            tmp_IDXxse{tmp_count}=['1/2*(-',num2str(tmp_Tlobe),'+',num2str(sym_x),')*2-(',tmp_IDXxse{tmp_count-4},')'];
            tmp_IDXyse{tmp_count}=tmp_IDXyse{tmp_count-4};
            tmp_count=tmp_count+1;
        end
        sym_x=sym_x+1;
    end
end
%point 2*Tlobe+3~4*Tlobe+3
for tmp_loop=(2*tmp_Tlobe+3):(4*tmp_Tlobe+3)
    tmp_IDXx{tmp_loop}=['-(',tmp_IDXx{2*(2*tmp_Tlobe+2)-tmp_loop},')'];
    tmp_IDXy{tmp_loop}=tmp_IDXy{2*(2*tmp_Tlobe+2)-tmp_loop};
    tmp_IDXxw{tmp_loop}=['-(',tmp_IDXxw{2*(2*tmp_Tlobe+2)-tmp_loop},')'];
    tmp_IDXyw{tmp_loop}=tmp_IDXyw{2*(2*tmp_Tlobe+2)-tmp_loop};
    tmp_IDXxse{tmp_loop}=['-(',tmp_IDXxse{2*(2*tmp_Tlobe+2)-tmp_loop},')'];
    tmp_IDXyse{tmp_loop}=tmp_IDXyse{2*(2*tmp_Tlobe+2)-tmp_loop};
end
%points 4*Tlobe+4~8*Tlobe+6
for tmp_loop=(4*tmp_Tlobe+4):(8*tmp_Tlobe+6)
    tmp_IDXx{tmp_loop}=tmp_IDXx{2*(4*tmp_Tlobe+3)-tmp_loop+1};
    tmp_IDXy{tmp_loop}=['-(',tmp_IDXy{2*(4*tmp_Tlobe+3)-tmp_loop+1},')'];
    tmp_IDXxw{tmp_loop}=tmp_IDXxw{2*(4*tmp_Tlobe+3)-tmp_loop+1};
    tmp_IDXyw{tmp_loop}=['-(',tmp_IDXyw{2*(4*tmp_Tlobe+3)-tmp_loop+1},')'];
    tmp_IDXxse{tmp_loop}=tmp_IDXxse{2*(4*tmp_Tlobe+3)-tmp_loop+1};
    tmp_IDXyse{tmp_loop}=['-(',tmp_IDXyse{2*(4*tmp_Tlobe+3)-tmp_loop+1},')'];
end
tmp_IDXx{tmp_Npts+1}=tmp_IDXx{1};
tmp_IDXy{tmp_Npts+1}=tmp_IDXy{1};
tmp_IDXxw{tmp_Npts+1}=tmp_IDXxw{1};
tmp_IDXyw{tmp_Npts+1}=tmp_IDXyw{1};
tmp_IDXxse{tmp_Npts+1}=tmp_IDXxse{1};
tmp_IDXyse{tmp_Npts+1}=tmp_IDXyse{1};

tmp_wallString=('''');
tmp_chloString=('''');
tmp_chliString=('''');
for tmp_loop=1:tmp_Npts
    tmp_wallString=[tmp_wallString,'(',tmp_IDXx{tmp_loop},')*scale_x'' ''(',tmp_IDXy{tmp_loop},')*scale_y'''];
    tmp_chloString=[tmp_chloString,'(',tmp_IDXxw{tmp_loop},')*scale_x'' ''(',tmp_IDXyw{tmp_loop},')*scale_y'''];
    tmp_chliString=[tmp_chliString,'(',tmp_IDXxse{tmp_loop},')*scale_x'' ''(',tmp_IDXyse{tmp_loop},')*scale_y'''];
    if tmp_loop~=tmp_Npts
        tmp_wallString=[tmp_wallString,char('; ''')];
        tmp_chloString=[tmp_chloString,char('; ''')];
        tmp_chliString=[tmp_chliString,char('; ''')];
    end
end
