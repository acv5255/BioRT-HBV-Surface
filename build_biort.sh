mkdir build && cd build
cmake ..
cmake --build . -- -j8
cp src/biort ..
rm -r build
