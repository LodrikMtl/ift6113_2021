
## Prerequisites
For ubuntu:
```bash
sudo apt-get install git build-essential cmake libx11-dev mesa-common-dev libgl1-mesa-dev libglu1-mesa-dev libxrandr-dev libxi-dev libxmu-dev libblas-dev libxinerama-dev libxcursor-dev
```
Other OS: see [libigl tutorial](https://libigl.github.io/tutorial/) and [libigl github issues](https://github.com/libigl/libigl/issues)


## Build
```c
mkdir build
cd build
cmake ..
make
./example
```
