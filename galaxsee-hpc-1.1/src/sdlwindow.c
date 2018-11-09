#ifdef HAS_SDL

#include "sdlwindow.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>

#include <SDL/SDL.h>
#include <SDL/SDL_opengl.h>
#include "quaternion.h"
#ifdef _REENTRANT
#ifdef _USE_PTHREADS
#include <pthread.h>
#endif
#else
#undef _USE_PTHREADS
#endif

#ifdef _USE_PTHREADS
pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;  
#endif

void sdlwindow_begin(sdlwindow * theWindow,int screenWidth,int screenHeight) {
    
    int error;
    error = SDL_Init(SDL_INIT_EVERYTHING);
    if(error) {
        /* ERROR HANDLING CODE HERE */
        /* (As if I would expect you to
         * care about error handling)   */
    }

    sdlwindow_init(theWindow,screenWidth,screenHeight);
}

void sdlwindow_init(sdlwindow * theWindow,int screenWidth,int screenHeight) {
    theWindow->anchored=0;
    Uint32 flags;
    int fullScreen     =    0; // Set to 1 for fullscreen
    double fieldOfView =   45; // In Degrees
    double nearClip    =    1; // How close must something be to not be drawn
    double farClip     =  200; // How far must something be to not be drawn

    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);

    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 16);
    SDL_GL_SetAttribute(SDL_GL_RED_SIZE  ,  8);
    SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE,  8);
    SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE ,  8);
    SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE,  8);

    flags = SDL_OPENGL;
    if(fullScreen) { flags |= SDL_FULLSCREEN; }
    theWindow->theContext =
        SDL_SetVideoMode(screenWidth, screenHeight, 0, flags);

    /* Initialization */
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);
    glEnable(GL_NORMALIZE);
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_COLOR_MATERIAL);
    //glShadeModel(GL_FLAT);

    /* Set up perspective */
    glViewport(0, 0, screenWidth, screenHeight);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(fieldOfView,
    (float)screenWidth/(float)screenHeight,
    nearClip, farClip);
    theWindow->rotator.a=-0.06;
    theWindow->rotator.b=0.8;
    theWindow->rotator.c=0.6;
    theWindow->rotator.d=0.3;

}

void sdlwindow_drawWorld(sdlwindow * theWindow) {
    double size = 1.0;
    double pos_x = 0.0;
    double pos_z = 0.0;
    if(theWindow->userDrawWorld!=NULL) {
        theWindow->userDrawWorld();
    } else {
        glColor3f(0.5,0.5,1.0);
        glBegin(GL_QUADS);
        glNormal3f(0.0,-1.0,0.0);
        glVertex3f(-size+pos_x,-size,-size+pos_z);
        glVertex3f(-size+pos_x,-size,size+pos_z);
        glVertex3f(size+pos_x,-size,size+pos_z);
        glVertex3f(size+pos_x,-size,-size+pos_z);
        glNormal3f(0.0,1.0,0.0);
        glVertex3f(-size+pos_x,size,-size+pos_z);
        glVertex3f(-size+pos_x,size,size+pos_z);
        glVertex3f(size+pos_x,size,size+pos_z);
        glVertex3f(size+pos_x,size,-size+pos_z);
        glEnd();
    }
}

void sdlwindow_render(sdlwindow * theWindow) {
    float param[4];
    fmatrix rotMatrix;

    /* Draw some stuff */
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // gluLookAt([Camera Location], [Camera Target], [Up Vector]);
    gluLookAt(
    0, 5, 1.5,
    0,  0, 0,
    0,  0, 1);


    // Diffusion color for LIGHT0
    // {red, green, blue, alpha}
    param[0] = 0.1; param[1] = 0.1; param[2] = 0.5; param[3] = 0.1;
    glLightfv(GL_LIGHT0, GL_AMBIENT, param);
    param[0] = 1.0; param[1] = 1.0; param[2] = 1.0; param[3] = 1.0;
    glLightfv(GL_LIGHT0, GL_DIFFUSE, param);

    // Location of LIGHT0
    // {x, y, z, [NOT USED]}
    param[0] = 1; param[1] = 1; param[2] = 1; param[3] = 0; 
    glLightfv(GL_LIGHT0, GL_POSITION, param);

    // Diffusion color for LIGHT1
    // {red, green, blue, alpha}
    param[0] = 0.1; param[1] = 0.1; param[2] = 0.5; param[3] = 0.1;
    glLightfv(GL_LIGHT1, GL_AMBIENT, param);
    param[0] = 1.0; param[1] = 1.0; param[2] = 1.0; param[3] = 1.0;
    glLightfv(GL_LIGHT1, GL_DIFFUSE, param);

    // Location of LIGHT1
    // {x, y, z, [NOT USED]}
    param[0] = -1; param[1] = -1; param[2] = -1; param[3] = 0; 
    glLightfv(GL_LIGHT1, GL_POSITION, param);

    glPushMatrix();
    QUAT_getMatrix(theWindow->rotator,&rotMatrix);
    glMultMatrixf(rotMatrix);
    // Put Geometry Drawing Code Here
#ifdef _USE_PTHREADS
    //if(!pthread_mutex_trylock(&mut)){
    pthread_mutex_lock(&mut);
#endif
    sdlwindow_drawWorld(theWindow);
    glPopMatrix();

    SDL_GL_SwapBuffers();
#ifdef _USE_PTHREADS
    pthread_mutex_unlock(&mut);
    //}
#endif
}

int sdlwindow_pollAndHandle(sdlwindow * theWindow) {
    SDL_Event event;

    int redraw=0;

    if(theWindow->userPollAndHandle!=NULL) {
        theWindow->userPollAndHandle();
    } else {
        SDL_PollEvent(&event);
        switch (event.type) {
            case SDL_MOUSEMOTION:
                //printf("Mouse moved by %d,%d to (%d,%d)\n", 
                //        event.motion.xrel, event.motion.yrel,
                //       event.motion.x, event.motion.y);
                if(theWindow->anchored) {
#ifdef _USE_PTHREADS
                    if(!pthread_mutex_trylock(&mut)) {
#endif
                    QUAT_mouseRotateSelf(&(theWindow->rotator),
                        1,3,0,0,event.motion.xrel,
                        event.motion.yrel,0.005);
                    redraw=1;
#ifdef _USE_PTHREADS
                    pthread_mutex_unlock(&mut);
                    }
#endif
                }
                break;
            case SDL_MOUSEBUTTONDOWN:
                //printf("Mouse button %d pressed at (%d,%d)\n",
                //       event.button.button, event.button.x, event.button.y);
                theWindow->anchored=1;
                redraw=1;
                break;
            case SDL_MOUSEBUTTONUP:
                //printf("Mouse button %d pressed at (%d,%d)\n",
                //       event.button.button, event.button.x, event.button.y);
                theWindow->anchored=0;
                redraw=1;
                break;
            case SDL_QUIT:
                exit(0);
                break;
            case SDL_KEYDOWN:
                if(event.key.keysym.sym==SDLK_q)
                    exit(0);
                break;
        }

    }
    return redraw;

}


void sdlwindow_end(sdlwindow * theWindow) {
    SDL_Quit();
}


#endif
