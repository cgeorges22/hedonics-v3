CXX=mpic++ -std=c++0x -Wno-unused-variable -Wno-unused-but-set-variable -Wno-maybe-uninitialized


queue_sim: simulation_helper.o hedonics-utilities-v2.10.o
	$(CXX) -o queue_sim simulation_helper.o hedonics-utilities-v2.10.o 

simulation_helper.o: simulation_helper.cpp
	$(CXX) -c  simulation_helper.cpp

hedonics-utilities-v2.10.o: hedonics-utilities-v2.10.cpp
	$(CXX) -c hedonics-utilities-v2.10.cpp

clean:
	rm queue_sim* simulation_helper.o hedonics-utilities-v2.10.o test.txt data* .data*

