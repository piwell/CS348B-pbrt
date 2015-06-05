#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_VOLUMES_RAINBOW_H
#define PBRT_VOLUMES_RAINBOW_H

// volumes/homogeneous.h*
#include "homogeneous.h"

// HomogenousDispersiveVolume Declarations
class RainbowVolume : public HomogeneousVolumeDensity {
public:
    // HomogeneousVolumeDensity Public Methods
    RainbowVolume(const Spectrum &sa, const Spectrum &ss, float gg,
        const Spectrum &emit, const BBox &e, const Transform &v2w):
        HomogeneousVolumeDensity(sa, ss, gg, emit, e, v2w) {
        
    }
    
    Spectrum waterdropReflection(const Spectrum& spectrum, const Vector &w, const Vector &wi);
        
private:

    float rainbowWavelength(const Vector &w, const Vector &wi);

};

RainbowVolume *CreateRainbowVolumeDensityRegion(const Transform &volume2world,
        const ParamSet &params);


#endif // PBRT_VOLUMES_RAINBOW_H
