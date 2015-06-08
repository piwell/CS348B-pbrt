
#include "stdafx.h"
#include "volumes/rainbow.h"
#include "paramset.h"


float LerpOrZero(float theta, float minTheta, float maxTheta,
                 float startWavelength, float endWavelength){

    if (theta < minTheta || maxTheta < theta){
        return 0;
    }

    const float thetaRange = maxTheta-minTheta;
    const float wavelengthRange = endWavelength-startWavelength;
    
    float wavelength = startWavelength + (theta-minTheta) * wavelengthRange / thetaRange;
    
    return wavelength;
}

float LerpTransfer(float x, float xMin, float xMax,
                 float y0, float y1){
    if (x < xMin){
        return y0;
    }
    if(xMax < x){
        return y1;
    }

    const float thetaRange = xMax-xMin;
    const float range = y1-y0;
    
    float y = y0 + (x-xMin) * range / thetaRange;
    
    return y;
}

//Wavelength and angle data from
//http://www.atoptics.co.uk/rainbows/sec.htm 
Spectrum RainbowVolume::rainbowReflection(const Spectrum& spectrum, 
            const Vector &w, const Vector &wi){
    //Phase angle
    float cosTheta = Dot(wi, -w);
    const float radToDeg = 57.2957;
    float theta = radToDeg * acosf(cosTheta);

    //Modified phase intensity for non-monocromatic light
    float I = PhaseMieHazy(wi, -w);
    float innerGlow = LerpTransfer(theta, 40.4, 40.45, 1.0, 0.9);
    I *= innerGlow;

    float rainbowI = 0.92f;
    float mistI = 0.08f;


    //Test primary rainbow
    float lambda = LerpOrZero(theta, 40.4, 42.3, 400.0, 700.0);
    if(!lambda) {
        //Test secondary rainbow
        lambda = LerpOrZero(theta, 51.0, 54.4, 700.0, 400.0);
        if(lambda){
            //Second rainbow has 42% intensity of primary
            rainbowI *= 0.42;
        }
    }
    if(!lambda) {
        return I * mistI * spectrum;
    }

    const Spectrum& rainbow = spectrum.filter(lambda);
    return I * (mistI*spectrum + rainbowI*rainbow);
}



RainbowVolume *CreateRainbowVolumeDensityRegion(const Transform &volume2world,
        const ParamSet &params){
    // Initialize common volume region parameters
    Spectrum sigma_a = params.FindOneSpectrum("sigma_a", 0.);
    Spectrum sigma_s = params.FindOneSpectrum("sigma_s", 0.);
    float g = params.FindOneFloat("g", 0.);
    Spectrum Le = params.FindOneSpectrum("Le", 0.);
    Point p0 = params.FindOnePoint("p0", Point(0,0,0));
    Point p1 = params.FindOnePoint("p1", Point(1,1,1));
    return new RainbowVolume(sigma_a, sigma_s, g, Le, BBox(p0, p1),
        volume2world);
}