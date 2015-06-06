
#include "stdafx.h"
#include "volumes/rainbow.h"
#include "paramset.h"

float RainbowVolume::rainbowWavelength(const Vector &w, const Vector &wi){
	float cosTheta = Dot(wi, -w);
    const float radToDeg = 57.2957;
    float theta = radToDeg * acosf(cosTheta);

	const float minTheta = 40.4, minWavelength = 400.0;
    const float maxTheta = 42.3, maxWavelength = 700.0;

	if (theta < minTheta || maxTheta < theta){
        return 0;
    }

	const float thetaRange = maxTheta-minTheta;
    const float wavelengthRange = maxWavelength-minWavelength;
    
    float wavelength = minWavelength + (theta-minTheta) * wavelengthRange / thetaRange;
    
    return wavelength;
}


Spectrum RainbowVolume::waterdropReflection(const Spectrum& spectrum, 
            const Vector &w, const Vector &wi){
    return spectrum.filter(rainbowWavelength(w, wi));
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