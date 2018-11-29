# CMake generated Testfile for 
# Source directory: /home/dennis/Documents/phd/courses/topics-in-computer-graphics/project/code/repo/pbrt-v3/src/ext/ptex/src/tests
# Build directory: /home/dennis/Documents/phd/courses/topics-in-computer-graphics/project/code/repo/pbrt-v3/build/src/ext/ptex/src/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(wtest "/home/dennis/Documents/phd/courses/topics-in-computer-graphics/project/code/repo/pbrt-v3/build/src/ext/ptex/src/tests/wtest")
add_test(rtest "/usr/bin/cmake" "-DOUT=/home/dennis/Documents/phd/courses/topics-in-computer-graphics/project/code/repo/pbrt-v3/build/src/ext/ptex/src/tests/rtest.out" "-DDATA=/home/dennis/Documents/phd/courses/topics-in-computer-graphics/project/code/repo/pbrt-v3/src/ext/ptex/src/tests/rtestok.dat" "-DCMD=/home/dennis/Documents/phd/courses/topics-in-computer-graphics/project/code/repo/pbrt-v3/build/src/ext/ptex/src/tests/rtest" "-P" "/home/dennis/Documents/phd/courses/topics-in-computer-graphics/project/code/repo/pbrt-v3/src/ext/ptex/src/tests/compare_test.cmake")
add_test(ftest "/usr/bin/cmake" "-DOUT=/home/dennis/Documents/phd/courses/topics-in-computer-graphics/project/code/repo/pbrt-v3/build/src/ext/ptex/src/tests/ftest.out" "-DDATA=/home/dennis/Documents/phd/courses/topics-in-computer-graphics/project/code/repo/pbrt-v3/src/ext/ptex/src/tests/ftestok.dat" "-DCMD=/home/dennis/Documents/phd/courses/topics-in-computer-graphics/project/code/repo/pbrt-v3/build/src/ext/ptex/src/tests/ftest" "-P" "/home/dennis/Documents/phd/courses/topics-in-computer-graphics/project/code/repo/pbrt-v3/src/ext/ptex/src/tests/compare_test.cmake")
add_test(halftest "/home/dennis/Documents/phd/courses/topics-in-computer-graphics/project/code/repo/pbrt-v3/build/src/ext/ptex/src/tests/halftest")
