#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CAMERAS_REALISTIC_H
#define PBRT_CAMERAS_REALISTIC_H

#include "pbrt.h"
#include "camera.h"
#include "film.h"
#include "film/image.h"

struct lensElement
{  
   public:
   float radius;
   float thickness;
   float zDist;
   float refraction;
   float aperture;
};

// Example representation of an autofocus zone.
class AfZone {
	public:
	  // from the data file
	  float left, right;
	  float top, bottom;
	  int xres,yres;
};
class RealisticCamera : public Camera {
public:
   RealisticCamera(const AnimatedTransform &cam2world,
      float hither, float yon, float sopen,
      float sclose, float filmdistance, float aperture_diameter,
      const string &specfile,
	  const string &autofocusfile,
      float filmdiag,
	  Film *film);
   ~RealisticCamera();
   float GenerateRay(const CameraSample &sample, Ray *) const;
   void  AutoFocus(Renderer * renderer, const Scene * scene, Sample * origSample);
   void  ParseAfZones(const string& filename);
   float computeF(Renderer * renderer, const Scene * scene, Sample * origSample,
                  AfZone* zone, float filmd);
private: 
   bool  autofocus;
   vector<AfZone> afZones;
   float ShutterOpen;
   float ShutterClose;
   Film * film;
   float filmdist;
   float filmDiag;
   vector<lensElement> lensElements;
};

RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film);

#endif
