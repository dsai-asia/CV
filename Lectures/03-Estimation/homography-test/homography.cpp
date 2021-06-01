/*
 * OpenCV homography estimation demo
 *
 * Simple demo program to load two images, get interest points, get correspondences
 * between interest points using SSD, get a homography between the images using RANSAC,
 * and display the results.
 *
 */

#include <iostream>
#include <opencv2/opencv.hpp>

using namespace cv;
using namespace std;

/* Names of the two images you want to process */

#define IMAGE1 "img1.png"
#define IMAGE2 "img2.png"

// Maximum number of interest points in each image

#define MAX_CORNERS 1000

// Distance ratio for best-to-worst distance for descriptor matches

#define DISTANCE_RATIO 4.0

// Maximum number of useful matches

#define MAX_MATCHES 500

/* Main: load images, get features and correspondences, estimate the homography, and show
 * the results. */

int main( int argc, char* argv[] ) {

    // Load images

    Mat matImage1 = imread(IMAGE1, IMREAD_GRAYSCALE);
    Mat matImage2 = imread(IMAGE2, IMREAD_GRAYSCALE);

    if (matImage1.rows == 0 || matImage1.cols == 0 || matImage2.rows == 0 || matImage2.cols == 0)
    {
        throw runtime_error("Could not load input images.");
    }

    // Display images

    namedWindow("Image 1", WINDOW_NORMAL);
    namedWindow("Image 2", WINDOW_NORMAL);

    imshow("Image 1", matImage1);
    imshow("Image 2", matImage2);

    // Get keypoints and keypoint descriptors for the two images

    vector<KeyPoint> vKeyPoints1;
    vector<KeyPoint> vKeyPoints2;

    Mat matDescriptors1;
    Mat matDescriptors2;

    Ptr<ORB> pOrb = ORB::create();
    pOrb->detectAndCompute(matImage1, Mat(), vKeyPoints1, matDescriptors1);
    pOrb->detectAndCompute(matImage2, Mat(), vKeyPoints2, matDescriptors2);

    // Do a brute force matching of the two point sets

    vector<DMatch> vMatches;
    BFMatcher matcher(NORM_L2, true);
    matcher.match(matDescriptors1, matDescriptors2, vMatches, Mat());

    // Sort matches by distance, remove weak matches, remove too many matches (snippet is from Learning OpenCV 3)

    sort(vMatches.begin(), vMatches.end());
    while (vMatches.front().distance * DISTANCE_RATIO < vMatches.back().distance)
    {
        vMatches.pop_back();
    }
    while (vMatches.size() > MAX_MATCHES)
    {
        vMatches.pop_back();
    }

    cout << "Got " << vMatches.size() << " matches." << endl;

    if (vMatches.size() < 4) {
        throw runtime_error("Insufficient matches to compute homography.");
    }

    vector<char> vCharMatchMask(vMatches.size(), 1);

    vector<Point2f> vPoints1;
    vector<Point2f> vPoints2;

    for (int i = 0; i < static_cast<int>(vMatches.size()); i++)
    {
        vPoints1.push_back(vKeyPoints1[vMatches[i].queryIdx].pt);
        vPoints2.push_back(vKeyPoints2[vMatches[i].trainIdx].pt);
    }
    findHomography(vPoints1, vPoints2, RANSAC, 4, vCharMatchMask);

    Mat matImageResult;
    drawMatches(matImage1, vKeyPoints1, matImage2, vKeyPoints2, vMatches, matImageResult, Scalar::all(-1),
                    Scalar::all(-1), vCharMatchMask, DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS);

    imshow("Matching result", matImageResult);

    waitKey(0);

    return 0;
}
