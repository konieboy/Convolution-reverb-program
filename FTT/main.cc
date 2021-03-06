// Sources:
// https://gist.github.com/tkaczenko/21ced83c69e30cfbc82b
// https://www.ee.iitb.ac.in/uma/~daplab/resources/wav_read_write.cpp
// Convolution algorithms from Leonard Manzara examples
// Tutorial examples

// How to profile
// - gcc -c filename
// double t1 = clock()
// Double finalTIme = clock() -1

// Overlap: decompose input to smaller signals and then combine components
// Discard the imaginary part of the 

#include <fstream>
#include <iostream>
#include <cstdio>
#include <cmath>
#include "string.h"
#include <sstream>
using namespace std;


#define SIZE       8
#define PI         3.141592653589793
#define TWO_PI     (2.0 * PI)
#define SWAP(a,b)  tempr=(a);(a)=(b);(b)=tempr

// Wave file Object
// Useful for loading a file, not really for saving one haha
class WaveFile {
    public:
        char ChunkID[4];
        int ChunkSize;
        char format[2];
        char subChunk1ID[4];
        int subChunk1Size;
        int16_t audioFormat;
        int16_t numChannels;
        int sampleRate;	
        int byteRate;
        int16_t blockAlign;
        int16_t bitsPerSample;
        char subChunk2ID[4];
        int subChunk2Size;
        float* fileData;

        int size;

   // Constructor with default values for data members
   WaveFile() {}
 
    void printHeader ()
    {
        cout << "File Data:" << endl;
        
        cout << "File Type: " << ChunkID << endl;
        cout << "File Size: " << ChunkSize << endl;
        cout << "WAV Marker: " << format << endl;
        cout << "Format Name: " << subChunk1ID << endl;
        cout << "Format Length: " << subChunk1Size << endl;
        cout << "Format Type: " << audioFormat << endl;
        cout << "Number of Channels: " << numChannels << endl;
        cout << "Sample Rate: " << sampleRate << endl;
        cout << "Sample Rate * Bits/Sample * Channels / 8: " << byteRate << endl;
        cout << "Bits per Sample * Channels / 8.1: " << blockAlign << endl;
        cout << "Bits per Sample: " << bitsPerSample << endl;
    }
};
// Code from Class handout:
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
float* convolve(float x[], int N, float h[], int M, float y[], int P)
{
    int n = 0, m = 0, t = 0;
    /*  Make sure the output buffer is the right size: P = N + M - 1  */
    if (P != (N + M - 1)) {
    printf("Output signal vector is the wrong size\n");
    printf("It is %-d, but should be %-d\n", P, (N + M - 1));
    printf("Aborting convolution\n");
    exit;
  }

  /*  Clear the output buffer y[] to all zero values  */  
    for (n = 0; n < P; n++)
        y[n] = 0.0;

    // /*  Do the convolution  */
    // /*  Outer loop:  process each input value x[n] in turn  */
    for (n = 0, t = 0; n < N; n++, t++) 
    {
        /*  Inner loop:  process x[n] with each sample of h[]  */
        for (m = 0; m < M; m++)
        {
            y[n+m] += x[n] * h[m];
        }
        if (t == 100000)
        {
            t = 0;
            cout << "x[" << n <<"]: " << x[n] << endl;
            cout << "h[" << m <<"]: " << h[m] << endl;

        }
    } 
    return y;
}

// Fast Convolve
//  The four1 FFT from Numerical Recipes in C,
//  p. 507 - 508.
//  Note:  changed float data types to double.
//  nn must be a power of 2, and use +1 for
//  isign for an FFT, and -1 for the Inverse FFT.
//  The data is complex, so the array size must be
//  nn*2. This code assumes the array starts
//  at index 1, not 0, so subtract 1 when
//  calling the routine (see main() below).
float* four1(float data[], int nn, int isign)
{
    unsigned long n, mmax, m, j, istep, i;
    float wtemp, wr, wpr, wpi, wi, theta;
    float tempr, tempi;

    n = nn << 1;
    j = 1;

    for (i = 1; i < n; i += 2) 
    {
        if (j > i) 
        {
            SWAP(data[j], data[i]);
            SWAP(data[j+1], data[i+1]);
        }

        m = nn;
        while (m >= 2 && j > m) 
        {
            j -= m;
            m >>= 1;
        }
        j += m;
    }

    mmax = 2;
    while (n > mmax) 
    {
	    istep = mmax << 1;
        theta = isign * (6.28318530717959 / mmax);
        wtemp = sin(0.5 * theta);
        wpr = -2.0 * wtemp * wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
	    for (m = 1; m < mmax; m += 2) 
        {
            for (i = m; i <= n; i += istep) 
            {
                j = i + mmax;
                tempr = wr * data[j] - wi * data[j+1];
                tempi = wr * data[j+1] + wi * data[j];
                data[j] = data[i] - tempr;
                data[j+1] = data[i+1] - tempi;
                data[i] += tempr;
                data[i+1] += tempi;
            }
            wr = (wtemp = wr) * wpr - wi * wpi + wr;
            wi = wi * wpr + wtemp * wpi + wi;
        }
	    mmax = istep;
    }
    return data;
}



// Creates a sine tone with the specified harmonic number.
// The array will be filled with complex numbers, and the
// signal is real (the imaginary parts are set to 0).

void createComplexSine(double data[], int size, int harmonicNumber)
{
    int i, ii;
    
    for (i = 0, ii = 0; i < size; i++, ii += 2) {
	data[ii] = sin((double)harmonicNumber * (double)i * TWO_PI / (double)size);
	data[ii+1] = 0.0;
    }
}



// Creates a cosine tone with the specified harmonic number.
// The array will be filled with complex numbers, and the
// signal is real (the imaginary parts are set to 0).

void createComplexCosine(double data[], int size, int harmonicNumber)
{
    int i, ii;

    for (i = 0, ii = 0; i < size; i++, ii += 2) {
	data[ii] = cos((double)harmonicNumber * (double)i * TWO_PI / (double)size);
	data[ii+1] = 0.0;
    }
}



// Creates a sawtooth wave, where each harmonic has
// the amplitude of 1 / harmonic_number.
// The array will be filled with complex numbers, and the
// signal is real (the imaginary parts are set to 0)

void createComplexSawtooth(double data[], int size)
{
    int i, ii, j;

    //  Calculate waveform using additive synthesis
    for (i = 0, ii = 0; i < size; i++, ii += 2) {
	data[ii] = 0.0;
	data[ii+1] = 0.0;
	for (j = 1; j <= size/2; j++) {
	    data[ii] +=
		(cos((double)j * (double)i * TWO_PI / (double)size)) / (double)j;
	}
    }
}



// Display the real and imaginary parts
// the data contained in the array.

void displayComplex(double data[], int size)
{
    int i, ii;

    printf("\t\tReal part \tImaginary Part\n");

    for (i = 0, ii = 0; i < size; i++, ii += 2)
	printf("data[%-d]: \t%.6f \t%.6f\n", i, data[ii], data[ii+1]);

    printf("\n");
}



// Performs the DFT on the input data,
// which is assumed to be a real signal.
// That is, only data at even indices is
// used to calculate the spectrum.

void complexDFT(double x[], int N)
{
    int n, k, nn;
    double omega = TWO_PI / (double)N;
    double *a, *b;

    // Allocate temporary arrays
    a = (double *)calloc(N, sizeof(double));
    b = (double *)calloc(N, sizeof(double));

    // Perform the DFT
    for (k = 0; k < N; k++) {
	a[k] = b[k] = 0.0;
	for (n = 0, nn = 0; n < N; n++, nn += 2) {
	    a[k] += (x[nn] * cos(omega * n * k));
	    b[k] -= (x[nn] * sin(omega * n * k));
	}
    }

    // Pack result back into input data array
    for (n = 0, k = 0; n < N*2; n += 2, k++) {
	x[n] = a[k];
	x[n+1] = b[k];
    }

    // Free up memory used for arrays
    free(a);
    free(b);
}



// Takes the results from a DFT or FFT, and
// calculates and displays the amplitudes of
// the harmonics.

double* postProcessComplex(double x[], int N)
{
    int i, k, j;
    double *amplitude, *result;

    // Allocate temporary arrays
    amplitude = (double *)calloc(N, sizeof(double));
    result = (double *)calloc(N, sizeof(double));

    // Calculate amplitude
    for (k = 0, i = 0; k < N; k++, i += 2) {
	// Scale results by N
	double real = x[i] / (double)N;
	double imag = x[i+1] / (double) N;
	// Calculate amplitude
	amplitude[k] = sqrt(real * real + imag * imag);
    }

    // Combine amplitudes of positive and negative frequencies
    result[0] = amplitude[0];
    result[N/2] = amplitude[N/2];
    for (k = 1, j = N-1; k < N/2; k++, j--)
	result[k] = amplitude[k] + amplitude[j];


    // Print out final result
    // printf("Harmonic \tAmplitude\n");
    // printf("DC \t\t%.6f\n", result[0]);
    // for (k = 1; k <= N/2; k++)
	// printf("%-d \t\t%.6f\n", k, result[k]);
    // printf("\n");

    // Free up memory used for arrays


    return result;
}




void saveWaveFile(string fileName, int outputFileSize, float* outputData, WaveFile inputData)
{
	ofstream outFile( fileName, ios::out | ios::binary);

	cout << endl << "Writing to file: " << fileName << endl;

	/*  Calculate the total number of bytes for the data chunk  */
	int chunkSize = inputData.numChannels * outputFileSize * (inputData.bitsPerSample / 8);

	char* chunkID = "RIFF";
	outFile.write( chunkID, 4);
	cout << "Chunk ID: " << chunkID << endl;

	outFile.write( (char*) &chunkSize, 4);
	cout << "Chunk Size: " << chunkSize << endl;
	cout << "Chunk Size OLD: " << inputData.ChunkSize << endl;


	char* format = "WAVE";
	outFile.write( format, 4);
	cout << "Format: " << format << endl;

	outFile.write( inputData.subChunk1ID, 4);
	cout << "Subchunk1ID: " << inputData.subChunk1ID << endl;

	int subChunk1Size = 16; // Might be 18?
	outFile.write( (char*) &subChunk1Size, 4);
	cout << "Subchunk1Size: " << subChunk1Size << endl;

	outFile.write( (char*) &inputData.audioFormat, 2);
	cout << "Audio Format: " << inputData.audioFormat << endl;

	outFile.write( (char*) &inputData.numChannels, 2);
	cout << "Number of Channels: " << inputData.numChannels << endl;

	outFile.write( (char*) &inputData.sampleRate, 4);
	cout << "sampleRate: " << inputData.sampleRate << endl;

	outFile.write( (char*) &inputData.byteRate, 4);
	cout << "ByteRate: " << inputData.byteRate << endl;

	outFile.write( (char*) &inputData.blockAlign, 2);
	cout << "BlockAlign: " << inputData.blockAlign << endl;

	outFile.write( (char*) &inputData.bitsPerSample, 2);
	cout << "BitsPerSample: " << inputData.bitsPerSample << endl;

	outFile.write( inputData.subChunk2ID, 4);
	cout << "Subchunk2ID: " << inputData.subChunk2ID << endl;

	//Data size
	int dataSize = outputFileSize * 2;
	outFile.write( (char*)&dataSize, 4);
	cout << "Data size: " << dataSize << " bytes" << endl;

	int16_t val;
	for(int i = 0; i < outputFileSize; i++)
	{
		val = (int16_t)(outputData[i] * 32767); // 2^15 - 1
        
		outFile.write((char*)&val, 2);
	}
	outFile.close();
}

WaveFile readFileData (string fileName)
{

    ifstream inputfile(fileName, ios::in | ios::binary); 
    inputfile.seekg(ios::beg);

    WaveFile fileData;

	inputfile.read(fileData.ChunkID, 4);

	inputfile.read((char*) &fileData.ChunkSize, 4);

	inputfile.read(fileData.format, 4);
	inputfile.read(fileData.subChunk1ID, 4);
	inputfile.read((char*) &fileData.subChunk1Size, 4);

	inputfile.read((char*) &fileData.audioFormat, 2);
	inputfile.read((char*) &fileData.numChannels, 2);
	inputfile.read((char*) &fileData.sampleRate, 4);
	inputfile.read((char*) &fileData.byteRate, 4);
	inputfile.read((char*) &fileData.blockAlign, 2);
	inputfile.read((char*) &fileData.bitsPerSample, 2);

	if (fileData.subChunk1Size == 18) 
    {
        char* spaceFiller = new char[2];
        inputfile.read( spaceFiller, 2);
    }

	inputfile.read(fileData.subChunk2ID, 4);
	inputfile.read((char*) &fileData.subChunk2Size, 4);
	fileData.size = fileData.subChunk2Size / 2;

	int size = fileData.size;
	
    cout << "size!" << fileData.size;

	short *data = new short[size];
	
   // short *data = new short[fileData.size];
    
	for (int i = 0; i < fileData.size; i++) 
    {
	    inputfile.read((char *) &data[i], 2);
    }
	
	inputfile.close();
	
	double sample; // Lots of precission baby

	float* dataOut = new float[fileData.size];

	printf("Size: %d\n", fileData.size);
	for (int i = 0; i < fileData.size; i++) 
    {
		sample = data[i];
		dataOut[i] = (sample) / (32767);

        // clamp to (-1,1) range to prevent pop (this doesn't actually seem to work)
        if (dataOut[i] > 1.0)
        {
			dataOut[i] = 1.0;
            cout << "clamped high";
        }
		if (dataOut[i] < -1.0)
        {
			dataOut[i] = -1.0;
            cout << "clamped low";
        }

        //dataOut[i] = 0.8 * dataOut[i];// scale to reduce pop

	}
    fileData.fileData = dataOut;
	return fileData;
}

// Multiply the contents of two arrays filled with complex doubles together
float* complexArrayMultiply (float X[], float H[], int length)
{
    double* output = new double[2 * length]; // X + H + 1

    for (int i = 0 ; i < length; i++)
    {
        // Get real component
        output[i] = ((X[i] * H[i]) - (X[i+1] * H[i+1]));
        // Get imaginary component
        output[i+1] = ((X[i] * H[i+1]) + (X[i+1] * H[i]));      
    }
    
    float* outputFloat = new float[2 * length];

    for (int i = 0 ; i < length; i++)
    {
        outputFloat[i] = (float) output[i];
    }

    return outputFloat;
}
// Code from http://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
// Find the upper bounded upper 2^n
int upper_power_of_two(int v)
{
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v++;
    return v;

}

int main ( int argc, char *argv[] )
{

    if ( argc < 4 ) 
    {
        cout<<"\n\nBad Args! Usage: convolve -inputWaveFile -IRfile -outputfile\n";
    }
    
    // Init Data
    string inputName = argv[1];
    string inpulseName = argv[2];
    string outputName = argv[3];
    
    // Read input WAV file
    cout <<"\n\nReading input data:\n\n" << endl;
    WaveFile inputData;
    inputData = readFileData(inputName);
    inputData.printHeader();


    // Read impulse WAV file
    cout <<"\n--------------------------------------------" << endl;
    cout <<"\nReading inpulse data:\n\n" << endl;
	WaveFile impuseData;
    impuseData = readFileData(inpulseName);
    impuseData.printHeader();


    // output
    int outputFileSize = (impuseData.size + inputData.size) - 1;

    // Slow code
        //float* outputFileData = new float[outputFileSize];
        
        //cout << "Convolving inputs...\n";
        //outputFileData = convolve(inputData.fileData, inputData.size, impuseData.fileData, impuseData.size, outputFileData, outputFileSize);
    
        //cout << "Saving output to: " << outputName;

        //saveWaveFile(outputName, outputFileSize,outputFileData,inputData);

    // Do fast convolution

    int sizeInputData = upper_power_of_two(2 * inputData.size);
    int arraySize = 2 * sizeInputData;

    float* paddedInputData = new float [arraySize]; // Array size has to be double the Upper power of 2
    float* paddedImpulseData = new float [arraySize];

    // Pad arrays with zero values
    // They must be the same size or else you get MONKEY TROUBLE
    for(int i = 0; i < arraySize; i++) 
    {
		paddedInputData[i] = 0.0;			
		paddedImpulseData[i] = 0.0;
	}	

    // Fill arrays with REAL part data
    // only copy to every seccond index
    for (int i = 0 ; i <  inputData.size; i++) // Input data
    {
        paddedInputData[i*2] = inputData.fileData[i];
    }

    for (int i = 0 ; i <  impuseData.size; i++) // Impulse data
    {
        paddedImpulseData[i*2] = impuseData.fileData[i];
    }

    // Get the Imaginary part of the data
    float* X = four1 (paddedInputData - 1, sizeInputData, 1);
    float* H = four1 (paddedImpulseData - 1, sizeInputData, 1);

    // Multiply together to get Y[n]
    float* complexRes = complexArrayMultiply(X, H, sizeInputData);

    // Run Reverse FTT on combined array
    float* finalFTT = four1 ((float*)complexRes - 1, sizeInputData, -1);

    // Run post processing at end
    // Scale down operation

    // Note: Unroll this later
    for (int i = 0 ; i < sizeInputData; i++) 
    {
        finalFTT[i] /= (float)sizeInputData;
    }

    // for (int i = 0 ; i <  inputData.size && i < 100; i++) // Input data
    // {
    //     cout << "finalFTT[i]: " << finalFTT[i] <<endl;
    // }

    // Save result to .wave    
    saveWaveFile(outputName, outputFileSize, finalFTT, inputData);

    return 0;
}

