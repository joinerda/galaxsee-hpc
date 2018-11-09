#ifndef SDLWINDOW
#define SDLWINDOW

#include "quaternion.h"
#include <SDL/SDL.h>
#include <SDL/SDL_opengl.h>


typedef struct {
    SDL_Surface * theContext;
    void (*userDrawWorld)();
    void (*userPollAndHandle)();
    Quaternion rotator;
    int anchored;
} sdlwindow;

void sdlwindow_begin(sdlwindow * thewindow,int screenWidth,int screenHeight);
void sdlwindow_init(sdlwindow * thewindow,int screenWidth,int screenHeight);
void sdlwindow_render(sdlwindow * thewindow);
void sdlwindow_drawWorld(sdlwindow * thewindow);
int sdlwindow_pollAndHandle(sdlwindow * thewindow);
void sdlwindow_end(sdlwindow * thewindow);

#endif
