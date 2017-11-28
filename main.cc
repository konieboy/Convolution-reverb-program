// Sources:
// https://gist.github.com/tkaczenko/21ced83c69e30cfbc82b
// https://www.ee.iitb.ac.in/uma/~daplab/resources/wav_read_write.cpp

// How to profile
// - gcc -c filename
// - Gcc 

#include <fstream>
#include <iostream>
#include <cstdio>
#include <cmath>
#include "string.h"

using namespace std;


/*****************************************************************************
*
*    Function:     convolve
*
*    Description:  Convolves two signals, producing an output signal.
*                  The convolution is done in the time domain using the
*                  "Input Side Algorithm" (see Smith, p. 112-115).
*
*    Parameters:   x[] is the signal to be convolved
*                  N is the number of samples in the vector x[]
*                  h[] is the impulse response, which is convolved with x[]
*                  M is the number of samples in the vector h[]
*                  y[] is the output signal, the result of the convolution
*                  P is the number of samples in the vector y[].  P must
*                       equal N + M - 1
*
*****************************************************************************/
void convolve(float x[], int N, float h[], int M, float y[], int P)
{
    int n, m;
    /*  Make sure the output buffer is the right size: P = N + M - 1  */
    if (P != (N + M - 1)) {
    printf("Output signal vector is the wrong size\n");
    printf("It is %-d, but should be %-d\n", P, (N + M - 1));
    printf("Aborting convolution\n");
    return;
  }
  /*  Clear the output buffer y[] to all zero values  */  
    for (n = 0; n < P; n++)
        y[n] = 0.0;
    /*  Do the convolution  */
    /*  Outer loop:  process each input value x[n] in turn  */
    for (n = 0; n < N; n++) 
    {
        /*  Inner loop:  process x[n] with each sample of h[]  */
        for (m = 0; m < M; m++)
        y[n+m] += x[n] * h[m];
    } 
}

// WAV Header data structure
struct wav_header
{
    char chunkID[4]; //"RIFF" = 0x46464952 (chunk descriptor)
    int chunkSize; //28 [+ sizeof(wExtraFormatBytes) + wExtraFormatBytes] + sum(sizeof(chunk.id) + sizeof(chunk.size) + chunk.size)
    char format[4]; //"WAVE" = 0x45564157

    char subchunk1ID[4]; //"fmt " = 0x20746D66 ! Must include space!
    int subchunk1Size; //16 [+ sizeof(wExtraFormatBytes) + wExtraFormatBytes]
    short int audioFormat;
    short int numChannels;
    int sampleRate;
    int byteRate;
    short int blockAlign;
    short int bitsPerSample;
    
 
    //[WORD wExtraFormatBytes;]
    //[Extra format bytes]
};

// Format of Chunks
struct chunk_t
{
    char ID[4]; //"data" = 0x61746164
    unsigned long size;  //Chunk data bytes
};

void printHeader (wav_header header)
{
     //Print WAV header
    printf("WAV File Header read:\n");

    
    cout << "File Type: " << header.chunkID << endl;
    cout << "File Size: " << header.chunkSize << endl;
    cout << "WAV Marker: " << header.format << endl;
    cout << "Format Name: " << header.subchunk1ID << endl;
    cout << "Format Length: " << header.subchunk1Size << endl;
    cout << "Format Type: " << header.audioFormat << endl;
    cout << "Number of Channels: " << header.numChannels << endl;
    cout << "Sample Rate: " << header.sampleRate << endl;
    cout << "Sample Rate * Bits/Sample * Channels / 8: " << header.byteRate << endl;
    cout << "Bits per Sample * Channels / 8.1: " << header.blockAlign << endl;
    cout << "Bits per Sample: " << header.bitsPerSample << endl;
}


int main ( int argc, char *argv[] )
{
    if ( argc < 2 ) 
        cout<<"\n\nBad Args. Usage: convolve inputfile IRfile outputfile\n";
    else 
    {
        FILE *inputFile = fopen(argv[1] , "rb");
        if (inputFile == NULL) {
            fprintf(stderr, "Can't open input file!\n");
            exit(1);
        }
        //Read WAV header
        wav_header header;
        fread(&header, sizeof(header), 1, inputFile);

        printHeader(header);
    }

}

