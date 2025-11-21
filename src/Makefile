CC = gcc
CFLAGS = -w -O3 -march=native
LDFLAGS = -lSDL2 -lm
TARGET = fluid_sim
SOURCE = sim.c

# --- Targets ---

# The default target.
all: $(TARGET)

# Rule to compile the executable
$(TARGET): $(SOURCE)
	$(CC) $(SOURCE) -o $(TARGET) $(CFLAGS) $(LDFLAGS)

# Rule to run the simulation after compile
run: $(TARGET)
	./$(TARGET)

# Rule to clean up the files
clean:
	rm -f $(TARGET)