# CMake generated Testfile for 
# Source directory: /home/dennis/Documents/phd/courses/topics-in-computer-graphics/project/code/repo/pbrt-v3
# Build directory: /home/dennis/Documents/phd/courses/topics-in-computer-graphics/project/code/repo/pbrt-v3/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(pbrt_unit_test "pbrt_test")
subdirs("src/ext/zlib")
subdirs("src/ext/openexr")
subdirs("src/ext/glog")
subdirs("src/ext/ptex")