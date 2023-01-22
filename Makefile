compile-examples: build-fft-example

warmup:
	-mkdir build\examples

build-fft-example: warmup
	$(CXX) examples/fft.cpp -g -o ./build/examples/fft -O2 -std=c++20 -Ifrpoly