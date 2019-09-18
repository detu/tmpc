bench: GemmNoWrapper.cpp Main.cpp Makefile
	g++ -std=c++17 -g -O2 -DNDEBUG -I/opt/blasfeo/include -L/opt/blasfeo/lib GemmNoWrapper.cpp Main.cpp -lblasfeo -lbenchmark -lpthread -lopenblas -o bench
all: bench
