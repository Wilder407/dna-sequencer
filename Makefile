CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Iinclude

all: nw

nw: app/main.cpp src/Matrix.cpp include/Matrix.h
	$(CXX) $(CXXFLAGS) app/main.cpp src/Matrix.cpp -o nw

test_matrix: tests/test_matrix.cpp src/Matrix.cpp include/Matrix.h
	$(CXX) $(CXXFLAGS) tests/test_matrix.cpp src/Matrix.cpp -o test_matrix

test: test_matrix
	./test_matrix

demo: nw
	./nw data/HBB_NM_000518.5_cds.fasta data/HBD_NM_000519.4_cds.fasta --shuffles 1000
	./nw data/HBB_NM_000518.5_cds.fasta data/HBB_sickle_c20A-T_cds.fasta --shuffles 1000
	./nw data/m_RNASequences.txt --free-end-gaps --shuffles 200

clean:
	rm -f nw test_matrix

.PHONY: all test demo clean
