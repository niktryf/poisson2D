### Compiler & flags
CC=gcc
CFLAGS=-c -fopenmp 

### Directories
SRC_DIR=src
OBJ_DIR=obj
HEAD_DIR=headers

### Executable name
EXEC=poisson2D

### Sources (in src directory)
SOURCES=$(addprefix $(SRC_DIR)/, \
	main.c io.c memory.c poisson2D.c \
	setup.c iterative.c)

### Objects (in obj directory)
OBJECTS=$(SOURCES:$(SRC_DIR)%.c=$(OBJ_DIR)%.o)

### Rules: #######################################

### General target (executable):
all: $(EXEC)

### How to make the executable:
$(EXEC): $(OBJECTS) 
	$(CC) $^ -o $@ -lm -fopenmp

### How to make every object:
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) $< -o $@

### How to clean up:
clean: 
	rm $(OBJECTS) $(EXEC)
