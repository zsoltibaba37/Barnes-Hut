
CC := g++
CCFLAGS := -Wfatal-errors -Wall -Wextra -std=c++17 -O3 -ffast-math -flto -fopenmp
LDFLAGS := -lm -lsfml-graphics -lsfml-window -lsfml-system

OUTPUT := Barnes-Hut
$(OUTPUT): $(wildcard *.cc)
	$(CC) $(CCFLAGS) $^ -o $@ $(LDFLAGS)

