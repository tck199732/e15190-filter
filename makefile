GCC := g++
CXX_FLAGS := -fPIC $(CXX_FLAGS) # path-independent code
CXX_FLAGS := -O2 $(CXX_FLAGS) # optimization
CXX_FLAGS := -std=c++17 $(CXX_FLAGS) # c++17
CXX_FLAGS := `root-config --cflags --libs` $(CXX_FLAGS) # root

main:
	$(GCC) main.cpp ./src/*.cpp -o main.exe  $(CXX_FLAGS) -I./include
