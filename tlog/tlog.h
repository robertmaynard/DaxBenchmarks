/*
 *  tlog.h
 *  parallelStreamline
 *
 *  Created by MiGi on 1/5/11.
 *  Copyright 2011 OSC. All rights reserved.
 *
 */

#ifndef TLOG_H
#define TLOG_H
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "pthread.h"
#include "cp_time.h"

class TLog{
    typedef unsigned long ulong;
    
    FILE *fp;
    int eventIDs;
    long long startTime;
    int threads;
    inline ulong getThreadID() {return (ulong)pthread_self();}
public:
    TLog(const char *filename_=NULL) : eventIDs(0), threads(0)
    {
#ifdef _PROFILE
        startTime = Timer::getTimeUS();
        
        char *filename;
        if (filename_==NULL) {
            filename = new char[1024];
            
            time_t rawtime;
            struct tm * info;
            time ( &rawtime );
            info = localtime ( &rawtime );
            sprintf ( filename, "tlog_%04d%02d%02d_%02d%02d%02d.txt", info->tm_year-100, info->tm_mon, info->tm_mday, info->tm_hour, info->tm_min, info->tm_sec );
            
        }else {
            filename = new char[strlen(filename_)+1];
            strcpy(filename, filename_);
        }
        
        fp = fopen(filename, "wt");
        if (fp)
            printf("Log filename: %s\n", filename);
        else {
            fprintf(stderr, "Cannot save to file %s!", filename);
        }

        delete[] filename;
#endif
    }
    ~TLog() 
    {
#ifdef _PROFILE
        fclose(fp);
#endif
    }
    
    
    void regThread(const char *name)
    {
#ifdef _PROFILE
        fprintf(fp, "T %lu %s\n", (ulong)pthread_self(), name);
#endif
    }
                
    int createEventID(const char *name, int r, int g, int b)
    {
#ifdef _PROFILE
        fprintf(fp, "D %d %d %d %d %s\n", eventIDs, r,g,b, name);
        return eventIDs++;
#else
        return 0;
#endif
    }
    
    ////////////
    inline void startEvent(int eventID)
    {
#ifdef _PROFILE
        fprintf(fp, "S %lld %lu %d\n", Timer::getTimeUS()-startTime, getThreadID(), eventID);
#endif
    }
    inline void startEvent(int eventID, const char *msg)
    {
#ifdef _PROFILE
        fprintf(fp, "S %lld %lu %d %s\n", Timer::getTimeUS()-startTime, getThreadID(), eventID, msg);
#endif
    }
    
    inline void endEvent(int eventID)
    {
#ifdef _PROFILE       
        fprintf(fp, "E %lld %lu %d\n", Timer::getTimeUS()-startTime, getThreadID(), eventID);
#endif
    }
    
    inline void logMessage(const char *msg)
    {
#ifdef _PROFILE
        fprintf(fp, "M %lld %lu %s\n", Timer::getTimeUS()-startTime, getThreadID(), msg);
#endif
    }
    
    
};

extern TLog *tlog;

#endif
