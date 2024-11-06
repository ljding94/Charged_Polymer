Compile=g++ -g -std=c++20 -O3 # improve run time e.g. 10s to 4s
charged_polymer: main.o charged_polymer.o
	$(Compile) -o charged_polymer main.o charged_polymer.o
main.o: main.cpp charged_polymer.h
	$(Compile) -o main.o -c main.cpp
charged_polymer.o: charged_polymer.cpp charged_polymer.h
	$(Compile) -o charged_polymer.o -c charged_polymer.cpp
clean:
	rm -f main.o charged_polymer.o nohup.out