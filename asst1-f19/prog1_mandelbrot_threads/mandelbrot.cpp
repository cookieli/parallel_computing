#include <stdio.h>
#include <pthread.h>

// Use this code to time your threads
#include "CycleTimer.h"


/*

  15418 Spring 2012 note: This code was modified from example code
  originally provided by Intel.  To comply with Intel's open source
  licensing agreement, their copyright is retained below.

  -----------------------------------------------------------------

  Copyright (c) 2010-2011, Intel Corporation
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    * Neither the name of Intel Corporation nor the names of its
      contributors may be used to endorse or promote products derived from
      this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
   IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
   TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
   PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
   OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


// Core computation of Mandelbrot set membershop
// Iterate complex number c to determine whether it diverges
static inline int mandel(float c_re, float c_im, int count)
{
    float z_re = c_re, z_im = c_im;
    int i;
    for (i = 0; i < count; ++i) {

        if (z_re * z_re + z_im * z_im > 4.f)
            break;

        float new_re = z_re*z_re - z_im*z_im;
        float new_im = 2.f * z_re * z_im;
        z_re = c_re + new_re;
        z_im = c_im + new_im;
    }

    return i;
}

//
// MandelbrotSerial --
//
// Compute an image visualizing the mandelbrot set.  The resulting
// array contains the number of iterations required before the complex
// number corresponding to a pixel could be rejected from the set.
//
// * x0, y0, x1, y1 describe the complex coordinates mapping
//   into the image viewport.
// * width, height describe the size of the output image
// * startRow, endRow describe how much of the image to compute
void mandelbrotSerial(
    float x0, float y0, float x1, float y1,
    int width, int height,
    int startRow, int endRow,
    int maxIterations,
    int output[])
{
    float dx = (x1 - x0) / width;
    float dy = (y1 - y0) / height;
    // printf("dx %f dy %f", dx, dy);
    for (int j = startRow; j < endRow; j++) {
        for (int i = 0; i < width; ++i) {
            float x = x0 + i * dx;
            float y = y0 + j * dy;

            int index = (j * width + i);
            output[index] = mandel(x, y, maxIterations);
        }
    }
    
}


// Struct for passing arguments to thread routine
typedef struct {
    float x0, x1;
    float y0, y1;
    unsigned int width;
    unsigned int height;
    int maxIterations;
    int* output;
    int threadId;
    int numThreads;

    float dx;
    float dy;
    int internal;
    int addition_count;
    int stage;
} WorkerArgs;



//
// workerThreadStart --
//
// Thread entrypoint.
void* workerThreadStart(void* threadArgs) {

    WorkerArgs* args = static_cast<WorkerArgs*>(threadArgs);

    int addition_per_thread = args->addition_count / args->numThreads;
    int additional_addition = args->addition_count % args->numThreads;
    int addition_start = args->height - args->addition_count;
    if(args->threadId < additional_addition){
        addition_start += (addition_per_thread+1) *args->threadId;
        addition_per_thread += 1;
    } else {
        addition_start  += (args->threadId - additional_addition)
                         * addition_per_thread
                         + (additional_addition)*(addition_per_thread + 1);
    }
    int write_count = args->stage + addition_per_thread;
    int cnt = 0;
    int addition_cnt = 0;
    for(int j = 0; j < write_count; j++){
        if(j != 0 && j % args->numThreads == 0){
            cnt++;
        }
        int start_row = (j % args->numThreads) * args->stage +args->threadId * args->internal + cnt;
        if(cnt >= args->internal) {
            start_row = addition_start + addition_cnt;
            addition_cnt++;
        }
        for(size_t i = 0; i < args->width; i++){
            float x = args->x0 + i* args->dx;
            float y = args->y0 + start_row * args->dy;

            int index = start_row * args->width + i;
            args->output[index] = mandel(x, y, args->maxIterations);
        }
    }
    //now we need test why it can not speed up linearly
    // double endTime = CycleTimer::currentSeconds();
    //  printf("this thread %d executes time [%.3f] ms\n", args->threadId, (endTime - startTime)*1000);
    return NULL;
}

//
// MandelbrotThread --
//
// Multi-threaded implementation of mandelbrot set image generation.
// Multi-threading performed via pthreads.
void mandelbrotThread(
    int numThreads,
    float x0, float y0, float x1, float y1,
    int width, int height,
    int maxIterations, int output[])
{
    const static int MAX_THREADS = 32;

    if (numThreads > MAX_THREADS)
    {
        fprintf(stderr, "Error: Max allowed threads is %d\n", MAX_THREADS);
        exit(1);
    }

    pthread_t workers[MAX_THREADS];
    WorkerArgs args[MAX_THREADS];
    float dx = (x1 - x0)/width;
    float dy = (y1 - y0)/height;
    int internal = height/(numThreads *numThreads);
    int addition_count = height%(numThreads * numThreads);
    int stage = internal * numThreads;

    for(int i = 0; i < numThreads; i++){
    // TODO: Set thread arguments here.
        args[i].dx = dx;
        args[i].dy = dy;
        args[i].internal = internal;
        args[i].addition_count = addition_count;
        args[i].stage = stage;
        args[i].threadId = i;
        args[i].x0 = x0;
        args[i].y0 = y0;
        args[i].x1 = x1;
        args[i].y1 = y1;
        args[i].width = width;
        args[i].height = height;
        args[i].output = output;
        args[i].maxIterations = maxIterations;
        args[i].numThreads = numThreads;
    }

    // Fire up the worker threads.  Note that numThreads-1 pthreads
    // are created and the main app thread is used as a worker as
    // well.

      for (int i=1; i<numThreads; i++)
        pthread_create(&workers[i], NULL, workerThreadStart, &args[i]);

    workerThreadStart(&args[0]);

    // wait for worker threads to complete
       for (int i=1; i<numThreads; i++)
          pthread_join(workers[i], NULL);
}
