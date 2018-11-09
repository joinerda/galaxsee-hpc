#ifndef READLINE_H
#define READLINE_H
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#define READLINE_MAX 1000

int stricmp_rl(char * word1, char * word2);
int readline(char * line, int maxline, FILE * fp);
int gettagline_rl(char * line,char * tag, char * rest);
void getint_rl(char * line,int * i);
void getdouble_rl(char * line, double * x);
void getword_rl(char * line, char * word);
#endif
