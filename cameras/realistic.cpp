// cameras/realistic.cpp*
#include "stdafx.h"
#include "cameras/realistic.h"
#include "paramset.h"
#include "sampler.h"
#include "montecarlo.h"
#include "filters/box.h"
#include "film/image.h"
#include "samplers/stratified.h"
#include "intersection.h"
#include "renderer.h"

#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <stdlib.h>

using namespace std;

RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film) {
	   // Extract common camera parameters from \use{ParamSet}
	   float hither = params.FindOneFloat("hither", -1);
	   float yon = params.FindOneFloat("yon", -1);
	   float shutteropen = params.FindOneFloat("shutteropen", -1);
	   float shutterclose = params.FindOneFloat("shutterclose", -1);

	   // Realistic camera-specific parameters
	   string specfile = params.FindOneString("specfile", "");
	   float filmdistance = params.FindOneFloat("filmdistance", 70.0); // about 70 mm default to film
	   float fstop = params.FindOneFloat("aperture_diameter", 1.0);
	   float filmdiag = params.FindOneFloat("filmdiag", 35.0);
	   string autofocusfile = params.FindOneString("af_zones", "");
	   assert(hither != -1 && yon != -1 && shutteropen != -1 &&
	      shutterclose != -1 && filmdistance!= -1);
	   if (specfile == "") {
	       Severe( "No lens spec file supplied!\n" );
	   }
	   return new RealisticCamera(cam2world, hither, yon,
	      shutteropen, shutterclose, filmdistance, fstop,
	      specfile, autofocusfile, filmdiag, film);
}

RealisticCamera::RealisticCamera(const AnimatedTransform &cam2world,
                                 float hither, float yon,
                                 float sopen, float sclose,
                                 float filmdistance, float aperture_diameter_,
                                 const string &specfile,
								 const string &autofocusfile,
                                 float filmdiag,
								 Film *f)
                                 : Camera(cam2world, sopen, sclose, f),
								   ShutterOpen(sopen),
								   ShutterClose(sclose),
								   film(f),
								   filmdist(filmdistance),
								   filmDiag(filmdiag)
{

	// YOUR CODE HERE -- build and store datastructures representing the given lens
	// and film placement.
	ifstream spec(specfile.c_str());
	if(!spec){
		fprintf(stderr, "Cannot open file %s\n", specfile.c_str());
      	exit (-1);
	}

	float zdist = 0;
	float thickness = 0;
	char line[512];
	while(!spec.eof()){
		spec.getline(line,512);
		 if (line[0] != '\0' && line[0] != '#' &&
         line[0] != ' ' && line[0] != '\t' && line[0] != '\n'){
		 	float r,z,n,a;
		 	sscanf(line, "%f %f %f %f\n",&r,&z,&n,&a);
		 	//float radius; float thickness; float zIntercept; float refraction; float aperture;
		 	a = (fabsf(r)>0) ? a:aperture_diameter_;

		 	lensElement elem;
		 	elem.radius 	= r;
		 	elem.aperture 	= a;
		 	elem.thickness 	= thickness;
		 	elem.refraction = n;
		 	elem.zDist 		= zdist;
		 	zdist -= z;
		 	thickness = z;
		 	//std::cout << zdist << std::endl;
		 	lensElements.push_back(elem);
         }
	}

	// If 'autofocusfile' is the empty string, then you should do
	// nothing in any subsequent call to AutoFocus()
	autofocus = false;

	if (autofocusfile.compare("") != 0)  {
		ParseAfZones(autofocusfile);
		autofocus = true;
	}
}

// parses the AF zone file
void RealisticCamera::ParseAfZones(const string& filename)
{
  ifstream specfile(filename.c_str());
   if (!specfile) {
      fprintf(stderr, "Cannot open file %s\n", filename.c_str());
      exit (-1);
   }

   char line[512];

   while (!specfile.eof()) {
      specfile.getline(line, 512);
      if (line[0] != '\0' && line[0] != '#' &&
         line[0] != ' ' && line[0] != '\t' && line[0] != '\n')
      {
		afZones.resize(afZones.size()+1);
		AfZone& zone = afZones[afZones.size()-1];
		sscanf(line, "%f %f %f %f\n", &zone.left, &zone.right, &zone.top, &zone.bottom);
      }
   }
	printf("Read in %zu AF zones from %s\n", afZones.size(), filename.c_str());
}

RealisticCamera::~RealisticCamera()
{

}

float RealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const
{
	//Raster to camera
	float rasX = sample.imageX-(film->xResolution/2);
	float rasY = sample.imageY-(film->yResolution/2);
	float scale = filmDiag/sqrt(pow(film->xResolution,2)+pow(film->yResolution,2));

	//Camera points
	float camX = -rasX*scale;
	float camY =  rasY*scale;

	lensElement firstLens 	= lensElements[lensElements.size()-1];
	float firstLensDist 	= firstLens.zDist;
	float fullFilmdist = firstLensDist-filmdist;

	if(firstLens.radius<0.f){
		float x = sqrt(pow(firstLens.radius,2)+pow(firstLens.aperture/2,2));
		firstLensDist += firstLens.radius+x;
	}

	//Lens points
	float lensU, lensV;
	ConcentricSampleDisk(sample.lensU,sample.lensV,&lensU,&lensV);
	float closestLensAperture = firstLens.aperture;
	lensU *= closestLensAperture/2;
	lensV *= closestLensAperture/2;

	Point pcam  = Point(camX,camY,fullFilmdist);
	Point plens = Point(lensU,lensV,firstLensDist);
	
	Ray lensRay = Ray(pcam,Normalize(plens-pcam),0.f,INFINITY);
	
	Vector v1 = Normalize(plens-pcam);
	Vector v2 = Vector(0.f,0.f,-1.f);
	float cosT = Dot(v1,v2);
	float w = M_PI*pow(firstLens.aperture/2,2)*(1/pow(filmdist,2))*pow(cosT,4);

	for(int i=lensElements.size()-1; i>=0; i--){
		lensElement elem = lensElements[i];

		float  R = elem.radius;
		Point  O = lensRay.o;
		Vector D = lensRay.d;

		if(R != 0.f){	
			Point  C = Point(0,0,elem.zDist-elem.radius);

			float a = Dot(D,D);
			float b = Dot(2*(O-C),D);
			float c = Dot(O-C,O-C)-R*R;
			float d = b*b-4*a*c;

			if(d<0.f){
				return 0.0f;
			}

			float t1 = (-b-sqrt(d))/(2*a);
			float t2 = (-b+sqrt(d))/(2*a);
			
			float t = 0.0f;
			if(t1>0.f && t2>0.f){
				t = (t1<t2) ? t1:t2;
			}else if(t2<0.f && t1>0.f){
				t = t1;
			}else if(t1<0.f && t2>0.f){
				t = t2;
			}	
				
			Point P = O+t*D;
			if(sqrt(P.x*P.x+P.y*P.y) > elem.aperture/2){
				return 0.0f;
			}
			

			Vector N = Normalize(P-C);
			if(R>0){
				N = -N;
			}

			float cosTh  = Dot(D,N);
			float n1 = elem.refraction;
			float n2 = (i!=0) ? lensElements[i-1].refraction : 1.0f;
			if(n2 == 0) n2 = 1;
			

			if(n1 != n2){
				float my = n1/n2;
				if(1-my*my*(1-cosTh*cosTh)<0){
					return 0.0f;
				}

				Vector T = my*D-(my*cosTh+sqrt(1-my*my*(1-cosTh*cosTh)))*N;
				lensRay.o = P;
				lensRay.d = Normalize(T);
			}
		}else{
			float t = (elem.zDist-lensRay.o.z)/(lensRay.d).z;
			Point P = lensRay.o+t*lensRay.d;

			if(sqrt(P.x*P.x+P.y*P.y) > elem.aperture/2){
				return 0.0f;
			}
		}
	}

	lensRay.o = lensRay.o;
    CameraToWorld(lensRay, ray);
    ray->d = Normalize(ray->d);
    ray->time = sample.time;

    return w;
}


//Help funxtions to computeF
int posInrgb(int i, int j, int width){
	return 3*(i+j*width);
}

float SML(float rgb[], int w, int h, int step, int N){
	float F=0;
	for(int x = step; x<w-step; x++){
		for(int y = step; y<h-step; y++){
			float rF = 	fabsf(2*rgb[posInrgb(x,y,w)  ]-rgb[posInrgb(x-step,y,w)  ]-rgb[posInrgb(x+step,y,w)  ])+
				 		fabsf(2*rgb[posInrgb(x,y,w)  ]-rgb[posInrgb(x,y-step,w)  ]-rgb[posInrgb(x,y+step,w)  ]);
			float gF = 	fabsf(2*rgb[posInrgb(x,y,w)+1]-rgb[posInrgb(x-step,y,w)+1]-rgb[posInrgb(x+step,y,w)+1])+
				 		fabsf(2*rgb[posInrgb(x,y,w)+1]-rgb[posInrgb(x,y-step,w)+1]-rgb[posInrgb(x,y+step,w)+1]);
			float bF = 	fabsf(2*rgb[posInrgb(x,y,w)+2]-rgb[posInrgb(x-step,y,w)+2]-rgb[posInrgb(x+step,y,w)+2])+
				 		fabsf(2*rgb[posInrgb(x,y,w)+2]-rgb[posInrgb(x,y-step,w)+2]-rgb[posInrgb(x,y+step,w)+2]);
			 F += (rF+gF+bF/3.f); 
		}
	}
	return F;
}

float RealisticCamera::computeF(Renderer * renderer, const Scene * scene, Sample * origSample, AfZone* zone, float filmd){
	float origFilmd = filmdist;
	filmdist = filmd;

	RNG rng;
	MemoryArena arena;
	Filter * filter = new BoxFilter(.5f,.5f);
	const float crop[] = {zone->left,zone->right,zone->top,zone->bottom};
	ImageFilm sensor(film->xResolution, film->yResolution, filter, crop,"foo.exr",false);
	int xstart,xend,ystart,yend;
	sensor.GetSampleExtent(&xstart,&xend,&ystart,&yend);

	StratifiedSampler sampler(xstart, xend, ystart, yend,
		                          16, 16, true, ShutterOpen, ShutterClose);

	// Allocate space for samples and intersections
	int maxSamples = sampler.MaximumSampleCount();
	Sample *samples = origSample->Duplicate(maxSamples);
	RayDifferential *rays = new RayDifferential[maxSamples];
	Spectrum *Ls = new Spectrum[maxSamples];
	Spectrum *Ts = new Spectrum[maxSamples];
	Intersection *isects = new Intersection[maxSamples];


	// Get samples from _Sampler_ and update image
	int sampleCount;
	while ((sampleCount = sampler.GetMoreSamples(samples, rng)) > 0) {
		// Generate camera rays and compute radiance along rays
		for (int i = 0; i < sampleCount; ++i) {
			// Find camera ray for _sample[i]_
			float rayWeight = this->GenerateRayDifferential(samples[i], &rays[i]);
			rays[i].ScaleDifferentials(1.f / sqrtf(sampler.samplesPerPixel));


			// Evaluate radiance along camera ray

			if (rayWeight > 0.f)
				Ls[i] = rayWeight * renderer->Li(scene, rays[i], &samples[i], rng,
												 arena, &isects[i], &Ts[i]);
			else {
				Ls[i] = 0.f;
				Ts[i] = 1.f;
			}

			// Issue warning if unexpected radiance value returned
			if (Ls[i].HasNaNs()) {
				Error("Not-a-number radiance value returned "
					  "for image sample.  Setting to black.");
				Ls[i] = Spectrum(0.f);
			}
			else if (Ls[i].y() < -1e-5) {
				Error("Negative luminance value, %f, returned"
					  "for image sample.  Setting to black.", Ls[i].y());
				Ls[i] = Spectrum(0.f);
			}
			else if (isinf(Ls[i].y())) {
				Error("Infinite luminance value returned"
					  "for image sample.  Setting to black.");
				Ls[i] = Spectrum(0.f);
			}

		}

		// Report sample results to _Sampler_, add contributions to image
		if (sampler.ReportResults(samples, rays, Ls, isects, sampleCount))
		{
			for (int i = 0; i < sampleCount; ++i)
			{
				sensor.AddSample(samples[i], Ls[i]);
			}
		}

		// Free _MemoryArena_ memory from computing image sample values
		arena.FreeAll();
	}

	float * rgb;
	int width;
	int height;
	sensor.WriteRGB(&rgb,&width,&height,1.f);
	// YOUR CODE HERE! The rbg contents of the image for this zone
	// are now stored in the array 'rgb'.  You can now do whatever
	// processing you wish
	float F = SML(rgb,width,height,2,5);
	//you own rgb  now so make sure to delete it:
	delete [] rgb;
	//if you want to see the output rendered from your sensor, uncomment this line (it will write a file called foo.exr)
	//sensor.WriteImage(1.f);


	delete[] samples;
	delete[] rays;
	delete[] Ls;
	delete[] Ts;
	delete[] isects;
	
	filmdist = origFilmd;
	return F;
}

void  RealisticCamera::AutoFocus(Renderer * renderer, const Scene * scene, Sample * origSample) {	
	if(!autofocus)
		return;

	float closestFocus = 0;
	for(size_t i=0; i<afZones.size(); i++){
		AfZone* zone = &afZones[i];

		float Fk1, Fk2, Fk3;
		float Fm1, Fm2, Fm3;
		Fm1 = Fm2 = Fm3 = 0;

		float deltaD = 5.f;
		float dm = std::max(filmdist-40.f,10.f);
		float bestd = 0;

		int M = 150.f/deltaD;
		int k = 2;

		Fk1 = computeF(renderer,scene,origSample,zone,dm-deltaD);
		Fk2 = computeF(renderer,scene,origSample,zone,dm);
		Fk3 = computeF(renderer,scene,origSample,zone,dm+deltaD);
			
		while(true){		
			if(Fk2>Fm2 && Fk2 > Fk3 && Fk2 > Fk1){
					Fm2   = Fk2;
					Fm1   = Fk1;
					Fm3   = Fk3;
					bestd = dm;
			}
			if(k>M){
				break;
			}

			dm = dm+deltaD;

			Fk1 = Fk2;
			Fk2 = Fk3;
			Fk3 = computeF(renderer,scene,origSample,zone,dm+deltaD);
			k++;
			if(Fm2 > 0 && (Fk1 + Fk2 + Fk3) < 10.f){
				break;
			}
		}

		float d = 	((log(Fm2)-log(Fm3))*(pow(bestd,2)-pow(bestd-deltaD,2)))-
				  	((log(Fm2)-log(Fm1))*(pow(bestd,2)-pow(bestd+deltaD,2)));
		d 		/=	2*deltaD*((log(Fm2)-log(Fm1))+(log(Fm2)-log(Fm3)));	

		closestFocus = std::max(d,closestFocus);
	}

	std::cout << "filmdist = " << closestFocus << std::endl;
	filmdist = closestFocus;
}
