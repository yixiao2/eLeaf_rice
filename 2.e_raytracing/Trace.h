/**************************************************************
*% eLeaf: 3D model of rice leaf photosynthesis
*% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
*% @author: Yi Xiao <yixiao20@outlook.com>
*% @version: 1.2.5
**************************************************************/

#ifndef TRACE_H
#define TRACE_H

extern void trace(Ray *ray,Object *obj);
extern double TriIsect(Ray *ray,Object *curobj,int i);
extern int JudgeInside(Ray *ray,Object *curobj,int i,double dis);

#endif

