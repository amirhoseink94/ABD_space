CC = g++
CFLAGS_NW = -g
CFLAGS = -Wall -g
OPGLFLAGS = -lGL -lGLU -lglut -lGLEW -lglfw -lX11 -lXxf86vm -lXrandr -lpthread -lXi -ldl -lXinerama -lXcursor
ARMAFLAG = -larmadillo
AVCFLAG = -lavcodec -lavformat -lavutil -lswscale
CVFLAGS = `pkg-config --cflags --libs opencv4`

res: main.o Particle.o Segment.o Face.o Body.o Engin.o Util.o
	$(CC) $(CFLAGS) -o res *.o $(OPGLFLAGS) $(ARMAFLAG) $(CVFLAGS)
main.o: main.cpp
	$(CC) $(CFLAGS) $(CVFLAGS) -c main.cpp
Util.o: src/Util.cpp
		$(CC) $(CFLAGS) -c src/Util.cpp
Engin.o: src/Engin.cpp
	$(CC) $(CFLAGS) -c src/Engin.cpp
Body.o: src/Body.cpp
	$(CC) $(CFLAGS) -c src/Body.cpp
Face.o: src/Face.cpp
	$(CC) $(CFLAGS) -c src/Face.cpp
Segment.o: src/Segment.cpp
	$(CC) $(CFLAGS) -c src/Segment.cpp
Particle.o: src/Particle.cpp
	$(CC) $(CFLAGS) -c src/Particle.cpp
clean:
	rm -f *.o
