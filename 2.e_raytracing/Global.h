/**************************************************************
*% eLeaf: 3D model of rice leaf photosynthesis
*% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
*% @author: Yi Xiao <yixiao20@outlook.com>
*% @version: 1.2.5
**************************************************************/

#ifndef GLOBAL_H
#define GLOBAL_H

extern void VecSub(Vec v1,Vec v2,Vec v3);	// vector from v1 to v2
extern void VecCross(Vec v1,Vec v2,Vec v3);	// v1 cross v2
extern void Normalize(Vec v);
extern void VecCopy(Vec v1,Vec v2);
extern double VecDot(Vec v1,Vec v2);
extern void IsectCopy(Isect *i1,Isect *i2);
extern void sum_p(void);

#endif
