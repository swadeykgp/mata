CC=g++

all: ass_scale ass_scale_heu ass_func_viz

ass_scale: ass_1_AIFA_scalability.cpp
	$(CC) -o ass_scale ass_1_AIFA_scalability.cpp

ass_scale_heu: ass_1_AIFA_scalability_only_heu.cpp
	$(CC) -o ass_scale_heu ass_1_AIFA_scalability_only_heu.cpp

ass_func_viz: ass_1_AIFA_functionality.cpp
	$(CC) -o ass_func_viz ass_1_AIFA_functionality.cpp

clean:
	rm -f  ass_func_viz ass_scale ass_scale_heu
