bassoon: main.cpp
	g++ -I /include -o bassoon main.cpp -O2 -larmadillo

clean:
	rm bassoon
