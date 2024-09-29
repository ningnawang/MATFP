# Computing Medial Axis Transform with Feature Preservation via Restricted Power Diagram :star2: SIGGRAPH Asia 2022

This project is licensed under the terms of the `MIT license`. Please refer to the following resources for more detail:

:cactus: [[**MATFP Project**](https://ningnawang.github.io/projects/2022_matfp/)]  :cactus: [[**MATFP Paper** (low res)](https://arxiv.org/abs/2210.13676)] :cactus: [[**MATFP Video** (YouTube)](https://youtu.be/0kP_EMtER-w?si=CyLhzKGUlTysoEUN)]  :cactus: [[**MATFP Video** (Bilibili)](https://www.bilibili.com/video/BV11hxqeqEcS/?share_source=copy_web&vd_source=085704da2cca04123412fb29bb28af85)]

![](./figures/teaser.png)

```
@article{2022MATFP,
  title={Computing Medial Axis Transform with Feature Preservation via Restricted Power Diagram},
  author={Wang, Ningna and Wang, Bin and Wang, Wenping and Guo, Xiaohu},
  journal={ACM Transactions on Graphics (Proceedings of SIGGRAPH Asia 2022)},
  volume={41},
  number={6},
  year={2022},
  publisher={ACM New York, NY, USA}
}
```

## Related Project MATTopo:
Our extended **SIGGRAPH Asia 2024** project **MATTopo: Topology-preserving Medial Axis Transform with Restricted Power Diagram** preserves the topological equivalence of the 3D medial mesh with respect to the input model.

:pineapple: [[**MATTopo Project**](https://ningnawang.github.io/projects/2024_mattopo/)] :pineapple: [[**MATTopo Paper** (low res)](https://arxiv.org/abs/2403.18761)] :pineapple: [[**MATTopo Code**](https://github.com/ningnawang/MATTopo)] :pineapple: [[**MATTopo Video** (YouTube)](https://www.youtube.com/watch?v=8AxJYVtU0SA)] :pineapple: [[**MATTopo Video** (Bilibili)](https://www.bilibili.com/video/BV1ZKxNeeEPF/?share_source=copy_web&vd_source=085704da2cca04123412fb29bb28af85)]


## Installation (MacOS)

We provide the commands for installing MATFP in MacOS:

- Clone the repository into your local machine:

```
git clone https://github.com/junanita/matfp.git
```

- Compile the code using cmake (first compilation may takes a while):

```
cd matfp
mkdir build && cd build
cmake ..
make -j4
```

You may need to install **gmp** and **mpfr** or **CGAL** before compiling the code. You can install them via [Homebrew](https://brew.sh/).

```
brew install gmp
brew install mpfr
brew install glfw
brew install cgal
```

Other external dependencies will be auto-downloaded during CMake compiling, including (verisons see `cmake/MATFPDownloadExternal.cmake`):
1. [cli11](https://github.com/CLIUtils/CLI11)
2. [fmt](https://github.com/fmtlib/fmt.git)
3. [geogram](https://github.com/alicevision/geogram.git)
4. [libigl](https://github.com/libigl/libigl.git)
5. [polyscope](https://github.com/nmwsharp/polyscope.git)
7. [spdlog](https://github.com/gabime/spdlog.git)


**NOTE**: During compilation, you may see error in geogram: `In grogram/src/bin/fpg/SymbolEnvironment.cpp "malloc.h" not found`. This is because "malloc.h" is deprecated, and can be simply fixed by removing line 11 "#include malloc.h".


## Usage 

The inputs of our software are triangle surface meshes in `.geogram` format with features (**sharp edges**, **concave edges** and **corners**) pre-detected and saved as [attributes](https://github.com/BrunoLevy/geogram/wiki/Mesh#attributes). We provide a light weight software to make this process a bit easier. Please note that `.geogram` format can also be viewed by [vorpaview](https://homepages.loria.fr/BLevy/GEOGRAM/vorpaview.html).

### Preprocess

The inputs of preprocess part `MATFP_PRE` are triangle surface meshes in `.off/.obj/.stl/.ply` format. Changing option `--concave` and `--convex` can manipulate the thresholds for detecting concave edges and sharp edges. And we assume that convex corners are adjacent to more than 2 sharp edges. We normalize all models in range `[0, 10]` by default. The models will be shown using [polyscope](https://github.com/nmwsharp/polyscope.git).


```
MATFP_PRE
Usage: ./MATFP_PRE [OPTIONS] input

Options:
  -h,--help                   Print this help message and exit
  --input TEXT REQUIRED       Input surface mesh INPUT in .off/.obj/stl/.ply format. (string, required)
  --sub INT                   Number of times for subdivision (unsigned int).
  --concave FLOAT             Threshold for detecting concave edges. Default 0.18, smaller more sensative (doouble)
  --convex FLOAT              Threshold for detecting sharp edges. Default 30, smaller more sensative (double)
  --save UINT                 Flag for saving preprocessed model (boolean).
```

For example:

```
$ ./MATFP_PRE ../models/69058.stl --sub=1 --concave=0.18 --convex=30 --save=1
```

will save prepocessed model as `input/mesh/mesh_69058.geogram`. 


![](https://github.com/junanita/matfp/blob/main/figures/matfp_pre.gif)


### Main usage

Once models are preprocessed and saved in `.geogram` format, we can run the main software. 

We offer two options to sample the surface seeds:
1. If the surface triangle is dense enough, we can sample centroids of triangles as surface seeds, and use option `--ds` to control the percentage of usage, we use `--ds=0.1` for most of models in [ABC Dataset](https://deep-geometry.github.io/abc-dataset/) since they are dense enough.
2. If the surface triangle is not dense, we can still sample random surface seeds by controlling [r-sample](https://www.cs.ucdavis.edu/~amenta/pubs/sm.pdf) as option `--r`. This option will override `--ds` and sample random surface seeds according to local feature size (LFS).

Please note that denser seeds may take longer time to process.


```
MATFP
Usage: ./MATFP [OPTIONS] input [downsample] [rsample]

Options:
  -h,--help                   Print this help message and exit
  --input TEXT REQUIRED       Input surface mesh INPUT in .geogram format. (string, required)
  --ds FLOAT                  Downsample percentage when input triangles are dense. Once set, we will not use rsample. Larger the denser. (double)
  --r FLOAT                   Control random sampling rate, r-sample from Nina Amenta (check power crust). Smaller the denser. (double)
  --cc FLOAT                  To specify pin points on concave edge with length cc_len_eps. Default=0.03. See Args.h for more detail. (double)
  ```


Example 1:

  ```
  $ ./MATFP ../input/mesh/mesh_00000068_767e4372b5f94a88a7a17d90_trimesh_009.geogram --ds=0.1
  ```
  
![](https://github.com/junanita/matfp/blob/main/figures/matfp_e1.gif)

Example 2:

For models with no sharp edges as MAT boundary, please setup thinning threshold (use `0.1` in gif).

```
$ ./MATFP ../input/mesh/mesh_bear.geogram --ds=0.1 
```

![](https://github.com/junanita/matfp/blob/main/figures/matfp_e2.gif)



### FAQ

1. ERROR during `make`: **'malloc.h' file not found**

It is common to see following error:

```
$PATH_TO_FOLDER/matfp/extern/geogram/src/bin/fpg/SymbolEnvironment.cpp:11:10: fatal error: 'malloc.h' file not found
```

The solution is easy, just comment line 11 under file `matfp/extern/geogram/src/bin/fpg/SymbolEnvironment.cpp`, such as:

```
# include <malloc.h>
```
Then rerun the `make -j4`.


## Copyright

All licenses in this repository are copyrighted by their respective authors. See `LICENSE` for details.
