CC = g++ 
CFLAGS = -O2 -I../c++ -mcx16 -march=native -fcilkplus -std=c++11 -Wall

index:	index.cpp
	$(CC) $(CFLAGS) index.cpp -o index

index2:	index2.cpp index.h
	$(CC) $(CFLAGS) index2.cpp -o index2

clean:
	rm -f index index2
