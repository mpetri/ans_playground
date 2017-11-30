
all:
	g++ -fcilkplus -O3 -DNDEBUG -std=c++14 -o ans_try.x ans_try.cpp
	g++ -fcilkplus -O0 -std=c++14 -o ans_try_dbg.x ans_try.cpp