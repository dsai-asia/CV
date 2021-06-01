
#include <cv.h>
#include <highgui.h>
#include <stdio.h>
#include <ctype.h>

/* Function prototypes */

IplImage *SkinDetect( IplImage *pImage );

/* Main function */

int main( int argc, char *argv[] ) {
    IplImage *pImage, *pImageSkin;
    char *pSzFilename;

    /* Get input image file name from command line */
    if ( argc != 2 ) {
        fprintf( stderr, "usage: %s <imagefile>\n", argv[0] );
        exit( -1 );
    }
    pSzFilename = strdup( argv[1] );

    pImage = cvLoadImage( pSzFilename, -1 );
    pImageSkin = SkinDetect( pImage );

    /* Display the image and wait for user to press a key */
    cvNamedWindow( "Testing OpenCV...", 1 );

    while( 1 ) {
        cvShowImage( "Testing OpenCV...", pImage );
        if ( tolower( cvWaitKey( 0 )) == 'q' ) break;
        cvShowImage( "Testing OpenCV...", pImageSkin );
        if ( tolower( cvWaitKey( 0 )) == 'q' ) break;
    }

    /* Clean up */
    cvDestroyWindow( "Testing OpenCV..." );
    cvReleaseImage( &pImage );
    cvReleaseImage( &pImageSkin );

    return 0;
}


/* Data structure for storing histogram information */

typedef struct {
    uchar cBins;
    double aHist[4096];
} tHistogram;


/* Read a histogram from a data file */

tHistogram *HistogramRead( char *pSzFilename ) {
    tHistogram *pHistogram;
    FILE *pFile;
    int cObj;

    /* Allocate memory */
    pHistogram = (tHistogram *)malloc( sizeof( tHistogram ));
    if ( !pHistogram ) return NULL;

    /* Open the file */
    pFile = fopen( pSzFilename, "rb" );
    if ( !pFile ) {
        fprintf( stderr, "%s %d: cannot open file %s\n",
                 __FILE__, __LINE__, pSzFilename );
        free( pHistogram );
        return NULL;
    }

    /* Read the histogram data */
    if (( cObj = fread( pHistogram, sizeof( tHistogram ), 1, pFile )) != 1 ) {
        fprintf( stderr, "%s %d: could not read sufficient data from %s\n",
                 __FILE__, __LINE__, pSzFilename );
        free( pHistogram );
        fclose( pFile );
        return NULL;
    }

    /* Finished */
    fclose( pFile );
    return pHistogram;
}


/* Free a histogram data structure */

void HistogramFree( tHistogram **ppHistogram ) {
    if ( *ppHistogram ) {
        free( *ppHistogram );
        *ppHistogram = NULL;
    }
}


/* Macros to get the max/min of 3 values */

#define MAX3(r,g,b) ((r)>(g)?((r)>(b)?(r):(b)):((g)>(b)?(g):(b)))
#define MIN3(r,g,b) ((r)<(g)?((r)<(b)?(r):(b)):((g)<(b)?(g):(b)))


/* Detect the skin pixels within an image */

IplImage *SkinDetect( IplImage *pImageIn ) {
    CvSize Size;
    int iRow, iCol;
    IplImage *pImageOut;
    uchar *aPixelIn, *aPixelOut;
    tHistogram *pHistSkin;
    tHistogram *pHistNonSkin;
    double fHistvalSkin, fHistvalNonSkin;

    /* Read skin and non-skin histograms from the data file */
    pHistSkin = HistogramRead( "poshist_hsv16bins.dat" );
    pHistNonSkin = HistogramRead( "neghist_hsv16bins.dat" );
    if ( pHistSkin == NULL || pHistNonSkin == NULL ) {
        fprintf( stderr, "%s %d: error: cannot read skin histogram data\n",
                 __FILE__, __LINE__ );
        exit( -1 );
    }

    /* Create an 8-bit image with one (gray) channel */
    Size.width = pImageIn->width;
    Size.height = pImageIn->height;
    pImageOut = cvCreateImage( Size, 8, 1 );

    /* Check that we have an RGB image */
    if ( pImageIn->nChannels != 3 ) {
        fprintf( stderr, "error: input image not RGB\n" );
        return pImageOut;
    }

    /* Classify the pixels */
    aPixelIn = (uchar *)pImageIn->imageData;
    aPixelOut = (uchar *)pImageOut->imageData;
    for ( iRow = 0; iRow < Size.height; iRow++ ) {
        for ( iCol = 0; iCol < Size.width; iCol++ ) {
            int R, B, G, F, I, X, H, S, V;

            /* Get RGB values -- OpenCV stores RGB images in BGR order!! */
            B = aPixelIn[ iRow * pImageIn->widthStep + iCol * 3 + 0 ];
            G = aPixelIn[ iRow * pImageIn->widthStep + iCol * 3 + 1 ];
            R = aPixelIn[ iRow * pImageIn->widthStep + iCol * 3 + 2 ];

            /* Convert RGB to HSV */
            X = MIN3( R, G, B );
            V = MAX3( R, G, B );
            if ( V == X ) {
                H = 0; S = 0;
            } else {
                S = (float)(V-X)/(float)V * 255.0;
                F = ( R==V ) ? (G-B) : (( G==V ) ? ( B-R ) : ( R-G ));
                I = ( R==V ) ? 0 : (( G==V ) ? 2 : 4 );
                H = ( I + (float)F/(float)(V-X) )/6.0*255.0;
                if ( H < 0 ) H += 255;
                if ( H < 0 || H > 255 || V < 0 || V > 255 ) {
                    fprintf( stderr, "%s %d: bad HS values: %d,%d\n",
                             __FILE__, __LINE__, H, S );
                    exit( -1 );
                }
            }

            /* Look up this H/S value in the skin/non-skin histogram tables */
            H = H >> 4;
            S = S >> 4;
            fHistvalSkin = pHistSkin->aHist[H*16+S];
            fHistvalNonSkin = pHistNonSkin->aHist[H*16+S];

            /* Set output pixel based on whether this is skin or not */
            if ( fHistvalSkin > fHistvalNonSkin ) {
                aPixelOut[ iRow * pImageOut->widthStep + iCol ] = V;
            } else {
                aPixelOut[ iRow * pImageOut->widthStep + iCol ] = 0;
            }
        }
    }

    /* Clean up */
    HistogramFree( &pHistSkin );
    HistogramFree( &pHistNonSkin );

    return pImageOut;
}

