CXX = g++
CXXFLAGS = -std=c++23

prototype:
	$(CXX) $(CXXFLAGS) -DHAHN_SERIES_DEMO prototype.cpp -o main
