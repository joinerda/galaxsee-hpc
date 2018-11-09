#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "readline.h"

int stricmp_rl(char * word1, char * word2) {
    int i;
    char temp1[READLINE_MAX];
    char temp2[READLINE_MAX];
    for(i=0;i<strlen(word1)+1;i++) temp1[i]=tolower(word1[i]);
    for(i=0;i<strlen(word2)+1;i++) temp2[i]=tolower(word2[i]);
    return strcmp(temp1,temp2);
}

int readline(char * line, int maxline, FILE * fp) {
    int line_read = 0;
    int length,i,loc;
    char temp[READLINE_MAX];
    char temp2[READLINE_MAX];
    while(!line_read) {
        if(fgets(temp,maxline,fp)==NULL) break;
        //strip comments
        strncpy(temp2,(const char *)temp,strcspn(temp,"#\n"));
        temp2[strcspn(temp,"#\n")]=0;
        strcpy(temp,temp2);
        length=strlen(temp);
        //strip trailing whitespace
        for(i=length-1;i>=0;i--) {
           if(temp[i]==' '||temp[i]=='\t'||temp[i]=='\n') {
               temp[i]='\0';
           } else break;
        }
        length=strlen(temp);
        loc=-1;
        //strip leading whitespace
        for(i=0;i<length;i++) {
            if(temp[i]==' '||temp[i]=='\t'||temp[i]=='\n') {
                loc=i;
            } else break;
        }
        if(loc>-1) {
            for(i=0;i<length-(loc);i++) {
                temp[i]=temp[i+loc+1];
            }
            temp[length-loc]='\0';
        }
        length=strlen(temp);
        if(length>0) line_read=1;
    }
    if(line_read) strcpy(line,temp);
    return line_read;  
} 

int gettagline_rl(char * line,char * tag, char * rest) {
    int loc,length,i;

    loc = strcspn(line," \t");
    strncpy(tag,line,loc);
    tag[loc]=0;
    length=strlen(line);
    for(i=0;i<length-(loc);i++) {
        rest[i]=line[i+loc+1];
    }
    rest[length-loc]='\0';
    return 1;
}
void getint_rl(char * line,int * i) {
    sscanf(line,"%d",i);
}
void getdouble_rl(char * line, double * x) {
    sscanf(line,"%lf",x);
}
void getword_rl(char * line, char * word) {
    sscanf(line,"%s",word);
}

