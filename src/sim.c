#include "sim.h"

int glb_windowWidth, glb_windowHeight,
    glb_cellSize, glb_gridWidth, glb_gridHeight,
    glb_globalDelay, glb_showGridLines, glb_generatedVelocity, glb_generatedDensity;

bool grid[MAX_GRID_H][MAX_GRID_W] = {false};  // true = liquid, false = empty

// Initialize SDL and create a window
SDL_Window *init_window() {
    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        fprintf(stderr, "SDL_Init Error: %s\n", SDL_GetError());
        return NULL;
    }

    SDL_Window *window = SDL_CreateWindow(
        WINDOW_NAME,
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
        glb_windowWidth, glb_windowHeight, SDL_WINDOW_SHOWN);

    if (!window) {
        fprintf(stderr, "SDL_CreateWindow Error: %s\n", SDL_GetError());
        SDL_Quit();
        return NULL;
    }

    return window;
}

//================================================================================================================================================================================================

// Simulation constants
#define DT 1.0 / 60.0
#define OVER_RELAXATION 1.9
#define DENSITY 1000.0
#define GRAVITY -9.81

// Grid dimensions
#define NUM_X 100
#define NUM_Y 100
#define NUM_CELLS (NUM_X * NUM_Y)
#define H (1.0 / NUM_Y)

// Field types
#define U_FIELD 0
#define V_FIELD 1
#define S_FIELD 2

typedef struct {
    float u[NUM_CELLS];
    float v[NUM_CELLS];
    float newU[NUM_CELLS];
    float newV[NUM_CELLS];
    float p[NUM_CELLS];
    float s[NUM_CELLS];  // smoke
    float m[NUM_CELLS];  // mass
    float newM[NUM_CELLS];
} Fluid;

typedef struct {
    Fluid fluid;
    int frameNr;
    bool showPressure;
    bool showSmoke;
    bool showVelocities;
    bool showStreamlines;
    bool paused;
    float obstacleX, obstacleY, obstacleRadius;
    bool showObstacle;
} Scene;

Scene scene;

bool SHOW_PRESSURE = false;
bool SHOW_SMOKE = true;

void initFluid() {
    memset(&scene.fluid, 0, sizeof(Fluid));

    // initialize all cells as fluids
    for (int i = 0; i < NUM_X; i++) {
        for (int j = 0; j < NUM_Y; j++) {
            scene.fluid.s[i * NUM_Y + j] = 1.0;

            if (SHOW_SMOKE)
                scene.fluid.m[i * NUM_Y + j] = 0;

            if (SHOW_PRESSURE)
                scene.fluid.m[i * NUM_Y + j] = 1;

            if (SHOW_SMOKE) {
                // add some initial smoke in the center
                if (i > NUM_X / 3 && i < 2 * NUM_X / 3 && j > NUM_Y / 3 && j < 2 * NUM_Y / 3) {
                    scene.fluid.m[i * NUM_Y + j] = 1.0;  // full smoke density
                }
            }

            // set boundaries
            if (i == 0 || i == NUM_X - 1 || j == 0) {
                scene.fluid.s[i * NUM_Y + j] = 0.0;
            }
        }
    }

    scene.showPressure = SHOW_PRESSURE;
    scene.showSmoke = SHOW_SMOKE;
    scene.paused = false;
    scene.frameNr = 0;
}

void initObstacle() {
    scene.obstacleX = NUM_X * H / 2.0f;  // Center position
    scene.obstacleY = NUM_Y * H / 2.0f;
    scene.obstacleRadius = 0.08f;
    scene.showObstacle = false;
}

void integrate(float dt) {
    for (int i = 1; i < NUM_X; i++) {
        for (int j = 1; j < NUM_Y - 1; j++) {
            if (scene.fluid.s[i * NUM_Y + j] != 0.0 &&
                scene.fluid.s[i * NUM_Y + j - 1] != 0.0) {
                scene.fluid.v[i * NUM_Y + j] += GRAVITY * dt;
            }
        }
    }
}

void solveIncompressibility(int numIters, float dt) {
    float cp = DENSITY * H / dt;

    for (int iter = 0; iter < numIters; iter++) {
        for (int i = 1; i < NUM_X - 1; i++) {
            for (int j = 1; j < NUM_Y - 1; j++) {
                if (scene.fluid.s[i * NUM_Y + j] == 0.0)
                    continue;

                float s = scene.fluid.s[i * NUM_Y + j];
                float sx0 = scene.fluid.s[(i - 1) * NUM_Y + j];
                float sx1 = scene.fluid.s[(i + 1) * NUM_Y + j];
                float sy0 = scene.fluid.s[i * NUM_Y + j - 1];
                float sy1 = scene.fluid.s[i * NUM_Y + j + 1];
                s = sx0 + sx1 + sy0 + sy1;
                if (s == 0.0)
                    continue;

                float div = scene.fluid.u[(i + 1) * NUM_Y + j] - scene.fluid.u[i * NUM_Y + j] +
                            scene.fluid.v[i * NUM_Y + j + 1] - scene.fluid.v[i * NUM_Y + j];

                float p = -div / s;
                p *= OVER_RELAXATION;
                scene.fluid.p[i * NUM_Y + j] += cp * p;

                scene.fluid.u[i * NUM_Y + j] -= sx0 * p;
                scene.fluid.u[(i + 1) * NUM_Y + j] += sx1 * p;
                scene.fluid.v[i * NUM_Y + j] -= sy0 * p;
                scene.fluid.v[i * NUM_Y + j + 1] += sy1 * p;
            }
        }
    }
}

void setObstacle(float x, float y, bool reset) {
    static bool firstCall = true;
    static float prevObstacleX = 0, prevObstacleY = 0;
    static float prevRadius = 0;

    // restore cells from previous obstacle position
    if (!firstCall && !reset) {
        float r = prevRadius;
        for (int i = 1; i < NUM_X - 2; i++) {
            for (int j = 1; j < NUM_Y - 2; j++) {
                float dx = (i + 0.5f) * H - prevObstacleX;
                float dy = (j + 0.5f) * H - prevObstacleY;

                // restore cells that are no longer in obstacle
                if (dx * dx + dy * dy < r * r) {
                    float new_dx = (i + 0.5f) * H - x;
                    float new_dy = (j + 0.5f) * H - y;
                    if (new_dx * new_dx + new_dy * new_dy >= scene.obstacleRadius * scene.obstacleRadius) {
                        scene.fluid.s[i * NUM_Y + j] = 1.0f;  // restore to fluid
                    }
                }
            }
        }
    }

    // mark new obstacle cels
    float r = scene.obstacleRadius;
    for (int i = 1; i < NUM_X - 2; i++) {
        for (int j = 1; j < NUM_Y - 2; j++) {
            float dx = (i + 0.5f) * H - x;
            float dy = (j + 0.5f) * H - y;

            if (dx * dx + dy * dy < r * r) {
                scene.fluid.s[i * NUM_Y + j] = 0.0f;  // mark as solid
                scene.fluid.u[i * NUM_Y + j] = 0.0f;
                scene.fluid.v[i * NUM_Y + j] = 0.0f;
                scene.fluid.m[i * NUM_Y + j] = 0.0f;  // clear smoke in obstacle
            }
        }
    }

    // store current position for next call
    prevObstacleX = x;
    prevObstacleY = y;
    prevRadius = r;
    firstCall = false;
}

void extrapolate() {
    for (int i = 0; i < NUM_X; i++) {
        scene.fluid.u[i * NUM_Y + 0] = scene.fluid.u[i * NUM_Y + 1];
        scene.fluid.u[i * NUM_Y + NUM_Y - 1] = scene.fluid.u[i * NUM_Y + NUM_Y - 2];
    }
    for (int j = 0; j < NUM_Y; j++) {
        scene.fluid.v[0 * NUM_Y + j] = scene.fluid.v[1 * NUM_Y + j];
        scene.fluid.v[(NUM_X - 1) * NUM_Y + j] = scene.fluid.v[(NUM_X - 2) * NUM_Y + j];
    }
}

float sampleField(float x, float y, int field) {
    x = fmaxf(fminf(x, NUM_X * H), H);
    y = fmaxf(fminf(y, NUM_Y * H), H);

    float dx = 0.0;
    float dy = 0.0;
    float *f;

    switch (field) {
    case U_FIELD:
        f = scene.fluid.u;
        dy = H * 0.5;
        break;
    case V_FIELD:
        f = scene.fluid.v;
        dx = H * 0.5;
        break;
    case S_FIELD:
        f = scene.fluid.m;
        dx = H * 0.5;
        dy = H * 0.5;
        break;
    }

    int x0 = (int)fminf(floorf((x - dx) / H), NUM_X - 1);
    float tx = ((x - dx) - x0 * H) / H;
    int x1 = fminf(x0 + 1, NUM_X - 1);

    int y0 = (int)fminf(floorf((y - dy) / H), NUM_Y - 1);
    float ty = ((y - dy) - y0 * H) / H;
    int y1 = fminf(y0 + 1, NUM_Y - 1);

    float sx = 1.0 - tx;
    float sy = 1.0 - ty;

    return sx * sy * f[x0 * NUM_Y + y0] +
           tx * sy * f[x1 * NUM_Y + y0] +
           tx * ty * f[x1 * NUM_Y + y1] +
           sx * ty * f[x0 * NUM_Y + y1];
}

void advectVel(float dt) {
    memcpy(scene.fluid.newU, scene.fluid.u, sizeof(scene.fluid.u));
    memcpy(scene.fluid.newV, scene.fluid.v, sizeof(scene.fluid.v));

    float h2 = 0.5 * H;

    for (int i = 1; i < NUM_X; i++) {
        for (int j = 1; j < NUM_Y; j++) {
            // u component
            if (scene.fluid.s[i * NUM_Y + j] != 0.0 &&
                scene.fluid.s[(i - 1) * NUM_Y + j] != 0.0 &&
                j < NUM_Y - 1) {
                float x = i * H;
                float y = j * H + h2;
                float u = scene.fluid.u[i * NUM_Y + j];
                float v = (scene.fluid.v[(i - 1) * NUM_Y + j] + scene.fluid.v[i * NUM_Y + j] +
                           scene.fluid.v[(i - 1) * NUM_Y + j + 1] + scene.fluid.v[i * NUM_Y + j + 1]) *
                          0.25;
                x = x - dt * u;
                y = y - dt * v;
                u = sampleField(x, y, U_FIELD);
                scene.fluid.newU[i * NUM_Y + j] = u;
            }
            // v component
            if (scene.fluid.s[i * NUM_Y + j] != 0.0 &&
                scene.fluid.s[i * NUM_Y + j - 1] != 0.0 &&
                i < NUM_X - 1) {
                float x = i * H + h2;
                float y = j * H;
                float u = (scene.fluid.u[i * NUM_Y + j - 1] + scene.fluid.u[i * NUM_Y + j] +
                           scene.fluid.u[(i + 1) * NUM_Y + j - 1] + scene.fluid.u[(i + 1) * NUM_Y + j]) *
                          0.25;
                float v = scene.fluid.v[i * NUM_Y + j];
                x = x - dt * u;
                y = y - dt * v;
                v = sampleField(x, y, V_FIELD);
                scene.fluid.newV[i * NUM_Y + j] = v;
            }
        }
    }

    memcpy(scene.fluid.u, scene.fluid.newU, sizeof(scene.fluid.u));
    memcpy(scene.fluid.v, scene.fluid.newV, sizeof(scene.fluid.v));
}

void advectSmoke(float dt) {
    memcpy(scene.fluid.newM, scene.fluid.m, sizeof(scene.fluid.m));

    float h2 = 0.5 * H;

    for (int i = 1; i < NUM_X - 1; i++) {
        for (int j = 1; j < NUM_Y - 1; j++) {
            if (scene.fluid.s[i * NUM_Y + j] != 0.0) {
                float u = (scene.fluid.u[i * NUM_Y + j] + scene.fluid.u[(i + 1) * NUM_Y + j]) * 0.5;
                float v = (scene.fluid.v[i * NUM_Y + j] + scene.fluid.v[i * NUM_Y + j + 1]) * 0.5;
                float x = i * H + h2 - dt * u;
                float y = j * H + h2 - dt * v;
                scene.fluid.newM[i * NUM_Y + j] = sampleField(x, y, S_FIELD);
            }
        }
    }
    memcpy(scene.fluid.m, scene.fluid.newM, sizeof(scene.fluid.m));
}

#define SPAWN_DENSITY 0.2f
#define SPAWN_VELOCITY 2.0f

void simulate() {
    if (scene.paused)
        return;

    if (scene.showObstacle) {
        setObstacle(scene.obstacleX, scene.obstacleY, false);
    }

    if (SHOW_PRESSURE) {
        for (int i = NUM_X / 4; i < 3 * NUM_X / 4; i++) {
            for (int j = NUM_Y / 4; j < 3 * NUM_Y / 4; j++) {
                scene.fluid.m[i + j] = SPAWN_DENSITY;   // full density
                scene.fluid.u[i + j] = SPAWN_VELOCITY;  // initial rightward velocity
            }
        }
    }
    /// SPAWNING SMOKE
    if (SHOW_SMOKE) {
        for (int i = NUM_X / 2 - 5; i < NUM_X / 2 + 5; i++) {
            for (int j = NUM_Y - 10; j < NUM_Y - 5; j++) {
                scene.fluid.m[i * NUM_Y + j] = 1.0;   // full smoke density
                scene.fluid.u[i * NUM_Y + j] = 0.0;   // no initial velocity
                scene.fluid.v[i * NUM_Y + j] = -2.0;  // up velocity
            }
        }

        for (int i = 5; i < 15; i++) {
            for (int j = NUM_Y - 10; j < NUM_Y - 5; j++) {
                scene.fluid.m[i * NUM_Y + j] = 1.0;
                scene.fluid.u[i * NUM_Y + j] = 1.0;
                scene.fluid.v[i * NUM_Y + j] = 0.0;
                // scene.fluid.m[i * NUM_Y + j] = 1.0;
                // scene.fluid.u[i * NUM_Y + j] = 0.0;
                // scene.fluid.v[i * NUM_Y + j] = -2.0;
            }
        }

        // for (int j = 50; j < 55; j++) {
        //     for (int i = NUM_X - 30; i < NUM_X - 25; i++) {
        //         scene.fluid.m[i * NUM_Y + j] = 1.0;
        //         scene.fluid.u[i * NUM_Y + j] = 1.0;
        //         scene.fluid.v[i * NUM_Y + j] = 0.0;
        //     }
        // }

        // for (int j = 50; j < 60; j++) {
        //     for (int i = NUM_X - 30; i < NUM_X - 25; i++) {
        //         scene.fluid.m[i * NUM_Y + j] = 1.0;
        //         scene.fluid.u[i * NUM_Y + j] = 1.0;
        //         scene.fluid.v[i * NUM_Y + j] = 0.0;
        //     }
        // }
    }

    integrate(DT);
    memset(scene.fluid.p, 0, sizeof(scene.fluid.p));
    solveIncompressibility(40, DT);
    extrapolate();
    advectVel(DT);
    advectSmoke(DT);
    // setObstacle(scene.obstacleX, scene.obstacleY, false);
    scene.frameNr++;
}

// Add smoke at a simulation coordinate
void addSmokeAt(float x, float y, float amount) {
    int i = (int)(x / H);
    int j = (int)(y / H);

    // Make sure we're within bounds
    if (i >= 0 && i < NUM_X && j >= 0 && j < NUM_Y) {
        // Only add smoke to fluid cells
        if (scene.fluid.s[i * NUM_Y + j] != 0.0f) {
            scene.fluid.m[i * NUM_Y + j] = fminf(scene.fluid.m[i * NUM_Y + j] + amount, 1.0f);
        }
    }
}

//================================================================================================================================================================================================

// Convert simulation value to color
SDL_Color getSciColor(float val, float minVal, float maxVal) {
    val = fmaxf(fminf(val, maxVal - 0.0001f), minVal);
    float d = maxVal - minVal;
    val = (d == 0.0f) ? 0.5f : (val - minVal) / d;

    float r, g, b;
    float m = 0.25f;
    int num = (int)(val / m);
    float s = (val - num * m) / m;

    switch (num) {
    case 0:
        r = 0.0f;
        g = s;
        b = 1.0f;
        break;
    case 1:
        r = 0.0f;
        g = 1.0f;
        b = 1.0f - s;
        break;
    case 2:
        r = s;
        g = 1.0f;
        b = 0.0f;
        break;
    case 3:
        r = 1.0f;
        g = 1.0f - s;
        b = 0.0f;
        break;
    default:
        r = g = b = 1.0f;
        break;
    }

    return (SDL_Color){r * 255, g * 255, b * 255, 255};
}

// Find min/max of pressure field (for visualization scaling)
void getPressureRange(float *minP, float *maxP) {
    *minP = *maxP = scene.fluid.p[0];
    for (int i = 0; i < NUM_CELLS; i++) {
        *minP = fminf(*minP, scene.fluid.p[i]);
        *maxP = fmaxf(*maxP, scene.fluid.p[i]);
    }
}

#define INNER_CIRCLE_R 0
#define INNER_CIRCLE_G 0
#define INNER_CIRCLE_B 0

#define OUTLINE_CIRCLE_R 255
#define OUTLINE_CIRCLE_G 255
#define OUTLINE_CIRCLE_B 255

#define CIRCLE_OUTLINE_SIZE 10

// Render the grid
void render_grid(SDL_Renderer *renderer) {
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderClear(renderer);

    float minP, maxP;
    getPressureRange(&minP, &maxP);

    for (int y = 0; y < NUM_Y; y++) {
        for (int x = 0; x < NUM_X; x++) {
            int idx = x * NUM_Y + y;

            SDL_Rect cell = {
                x * DEF_CELL_SIZE,
                y * DEF_CELL_SIZE,
                DEF_CELL_SIZE,
                DEF_CELL_SIZE};

            // Choose what to display based on visualization mode
            if (scene.showPressure) {
                float p = scene.fluid.p[idx];
                SDL_Color color = getSciColor(p, minP, maxP);
                if (scene.showSmoke) {
                    float s = scene.fluid.m[idx];
                    color.r = fmaxf(0, color.r - 255 * s);
                    color.g = fmaxf(0, color.g - 255 * s);
                    color.b = fmaxf(0, color.b - 255 * s);
                }
                SDL_SetRenderDrawColor(renderer, color.r, color.g, color.b, 255);
            } else if (scene.showSmoke) {
                float s = scene.fluid.m[idx];
                SDL_SetRenderDrawColor(renderer,
                                       150 + 105 * s,
                                       180 + 75 * s,
                                       255,
                                       255);
            } else if (scene.fluid.s[idx] == 0.0f) {
                // Solid boundary
                SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
            } else {
                // Empty fluid cell
                SDL_SetRenderDrawColor(renderer, 50, 50, 50, 255);
            }

            SDL_RenderFillRect(renderer, &cell);

            // Optional: Draw grid lines
            if (glb_showGridLines) {
                SDL_SetRenderDrawColor(renderer, 20, 20, 20, 255);
                SDL_RenderDrawRect(renderer, &cell);
            }
        }
    }

    if (scene.showObstacle) {
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);  // Gray obstacle
        int centerX = (int)(scene.obstacleX * DEF_CELL_SIZE / H);
        int centerY = (int)(scene.obstacleY * DEF_CELL_SIZE / H);
        int radius = (int)(scene.obstacleRadius * DEF_CELL_SIZE / H);
        int radiusSquared = radius * radius;
        // Draw filled circle
        for (int y = -radius; y <= radius; y++) {
            for (int x = -radius; x <= radius; x++) {
                int currentPoint = x * x + y * y;
                if (currentPoint < radius * radius) {
                    SDL_SetRenderDrawColor(renderer, INNER_CIRCLE_R, INNER_CIRCLE_G, INNER_CIRCLE_B, 255);
                    SDL_RenderDrawPoint(renderer, centerX + x, centerY + y);
                }
            }
        }
    }

    SDL_RenderPresent(renderer);
}

//================================================================================================================================================================================================

void get_new_setting(char *answer, int *setting) {
    *setting = atoi(answer);
}

int check_new_setting(int newSetting, int min, int max) {
    if (newSetting < min || newSetting > max) {
        printf("Wrong input, the value you enter must be between %d and %d!\n", min, max);
        printf("Please enter it again.\n");
        return 0;
    }

    return 1;
}

void process_new_setting(char *answer, const char *settingText, int *settingAddrs, int min, int max) {
    int passed = 0;
    while (!passed) {
        printf("%s", settingText);
        scanf("%s", answer);
        get_new_setting(answer, settingAddrs);
        passed = check_new_setting(*settingAddrs, min, max);
    }

    // printf("\n");
}

void print_settings() {
    printf("================================SETTINGS================================\n");
    printf("%s", SETTING_1_STRING);
    printf("%d\n", glb_windowHeight);

    printf("%s", SETTING_2_STRING);
    printf("%d\n", glb_windowWidth);

    printf("%s", SETTING_3_STRING);
    printf("%d\n", glb_cellSize);

    printf("%s", SETTING_4_STRING);
    printf("%d\n", glb_globalDelay);

    printf("%s", SETTING_5_STRING);
    printf("%d\n", glb_showGridLines);

    printf("========================================================================\n");
}

void settings_start() {

    print_settings();

    // printf("Window height : %d\n", DEF_WINDOW_H);
    // printf("Window width : %d\n", DEF_WINDOW_W);
    // printf("Grid cell size : %d\n", glb_cellSize);
    // printf("Delay (inverse frame rate) : %d\n", DEF_GLOBAL_DELAY);

    printf("\nWould you like to change these settings? [Y/N]\n");

    char *answer = malloc(100);
    scanf("%s", answer);
    if (answer[0] == 'Y' || answer[0] == 'y') {
        printf("Please introduce the new settings. \n\n");

        process_new_setting(answer, SETTING_1_STRING, &glb_windowWidth,
                            MIN_WINDOW_SIZE, MAX_WINDOW_SIZE);
        process_new_setting(answer, SETTING_2_STRING, &glb_windowHeight,
                            MIN_WINDOW_SIZE, MAX_WINDOW_SIZE);
        process_new_setting(answer, SETTING_3_STRING, &glb_cellSize,
                            MIN_CELL_SIZE, MAX_CELL_SIZE);
        process_new_setting(answer, SETTING_4_STRING, &glb_globalDelay,
                            MIN_GLOBAL_DELAY, MAX_GLOBAL_DELAY);
    }

    // printf("\n\nChoose what simulation to start ('P' for pressure, default is smoke)\n\n");
    // scanf("%s", answer);

    // if (answer[0] == 'P') {
    //     SHOW_PRESSURE = true;
    //     SHOW_SMOKE = false;
    // }

    // print_settings();

    // sleep(4);

    free(answer);
    return;
}

// void initialize_grid_dimensions() {
//     glb_gridWidth = glb_windowWidth / glb_cellSize;
//     glb_gridHeight = glb_windowWidth / glb_cellSize;
// }

void initialize_vars() {
    glb_windowWidth = DEF_WINDOW_W;
    glb_windowHeight = DEF_WINDOW_H;
    glb_cellSize = DEF_CELL_SIZE;
    glb_globalDelay = DEF_GLOBAL_DELAY;
}

int main(int argc, char *argv[]) {

    initialize_vars();

    if (argc > 1 && !strcmp(argv[1], "-s"))
        settings_start();
    // initialize_grid_dimensions();

    SDL_Window *window = init_window();
    if (!window)
        return 1;

    SDL_Renderer *renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    if (!renderer) {
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }

    initFluid();
    initObstacle();
    bool running = true;
    SDL_Event event;

    while (running) {
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                running = false;
            } else if (event.type == SDL_KEYDOWN) {
                // printf("Pressed key!\n");
                switch (event.key.keysym.sym) {
                case SDLK_o:
                    scene.showObstacle = !scene.showObstacle;
                    break;
                case SDLK_UP:
                    scene.obstacleY -= 0.05f;
                    break;
                case SDLK_DOWN:
                    scene.obstacleY += 0.05f;
                    break;
                case SDLK_LEFT:
                    scene.obstacleX -= 0.05f;
                    break;
                case SDLK_RIGHT:
                    scene.obstacleX += 0.05f;
                    break;
                }
            }
        }

        simulate();                  // Update fluid
        render_grid(renderer);       // Draw grid
        SDL_Delay(glb_globalDelay);  // Control simulation speed
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}