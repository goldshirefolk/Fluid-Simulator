#include <SDL2/SDL.h>
#include <stdbool.h>
#include <stdio.h>
#include <unistd.h>

#define WINDOW_NAME "Fluid Simulator"

/// DEF -> default
#define DEF_WINDOW_W 800
#define DEF_WINDOW_H 800
#define MIN_WINDOW_SIZE 200
#define MAX_WINDOW_SIZE 1600

#define DEF_CELL_SIZE 10
#define MIN_CELL_SIZE 2
#define MAX_CELL_SIZE 30
// #define DEF_GRID_W (DEF_WINDOW_W / DEF_CELL_SIZE)
// #define DEF_GRID_H (DEF_WINDOW_H / DEF_CELL_SIZE)
#define MAX_GRID_W 800
#define MAX_GRID_H 800
#define MIN_GLOBAL_DELAY 10
#define MAX_GLOBAL_DELAY 200

#define DEF_GLOBAL_DELAY 10

#define SETTING_1_STRING "Window height : "
#define SETTING_2_STRING "Window width : "
#define SETTING_3_STRING "Grid cell size : "
#define SETTING_4_STRING "Delay (inverse frame rate) : "
#define SETTING_5_STRING "Show grid lines : "