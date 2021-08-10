# FF
This is an implementation of the paper "Practical Foldover-Free Volumetric Mapping Constructio".

### Dependencies
* [OpenMesh](https://www.graphics.rwth-aachen.de/software/openmesh/download/)
* [OpenVolumeMesh](https://www.graphics.rwth-aachen.de/software/openvolumemesh/download/)
* [Eigen](http://eigen.tuxfamily.org/)
* [PARDISO](https://pardiso-project.org/)
* [Qt](http://download.qt.io/archive/qt/) optional for GUI

### Usage

```
2D case:
source is a planar mesh: PP2D_code.exe source.obj init.obj result.obj
sourc is a surface mesh: PP2D_code.exe source.obj init.obj result.obj

3D case:
PP_Code.exe source_mesh.ovm initial_mesh.ovm
```
Note that: tsource_mesh.ovm is source tetrahedral mesh; initial_mesh.ovm is initial volumetric map which contain foldovers;
The parameter `num_procs` should be set according to your CPU cores number.
