# Three Dimensional Modelling of Rice Leaf Photosynthesis
eLeaf is an end-to-end 3D modelling tool for rice leaf photosynthesis. Users just need to input a list of anatomical and biochemical parameters, run a single line of code in MATLAB, then eLeaf will automatically 1) reconstruct a representative 3D leaf anatomy, 2) simulate the light absorption profile with a raytracing algorithm, 3) simulate the reaction-diffusion system of CO2 inside the leaf. 

**[Reference]**
Xiao Y, Sloan J, Hepworth C, Fradera-Soler M, Mathers A, Thorley R, Baillie A, Jones H, Chang TG, Osborne CP, Chen XY, Sturrock C, Mooney S, Fleming AJ, Zhu X-G. Identifying leaf structure for improved rice performance via eLeaf, an integrated spatial model of photosynthesis. (2021) under review.

## Getting Started
In eLeaf model, 3D reconstruction was achieved through COMSOL Multiphysics and MATLAB. Ray tracing was coded in C program. Reaction-diffusion was simulated and solved by COMSOL Multiphysics. All modules were mastered by several MATLAB scripts.  

### - Prerequisites
Here we used a tower workstation (2x Intel Xeon E5 2650v4; 48 Cores; 128GB RAM) with CentOS 7 operating system.

1. COMSOL Multiphyscis:  
   - COMSOL version >= 5.3
   - Required modules of COMSOL:   
     - CFD module  
     - CAD import module  
     - LiveLink for MATLAB module  

2. MATLAB:  
   - Any compatible MATLAB, we used MATLAB 2015b.  

3. gcc:  
   - We used 4.8.5.  

## Versioning

This is eLeaf code version 1.2.5, first released on Jan. 27th 2021.

## Authors

* code by **Yi Xiao** - *yixiao20@outlook.com* - [https://github.com/xiaoyizz78](https://github.com/xiaoyizz78)

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.md](LICENSE.md) file for details

