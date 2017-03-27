#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
// Needed on MsWindows
#define NOMINMAX
#include <windows.h>
#endif // Win32 platform

#include <GLUT/glut.h>
// Download glut from: http://www.opengl.org/resources/libraries/glut/

#include "float.h"
#include "float2.h"
#include "float3.h"
#include "float4.h"
#include "float4x4.h"
#include <vector>

// MATERIALS

class Material
{
public:
    virtual float3 shade(
                         float3 normal,
                         float3 viewDir,
                         float3 lightDir,
                         float3 lightPowerDensity,
                         float3 position) {
        return float3(0,0,0);
    }
    float3 snoiseGrad(float3 r)
    {
        unsigned int x = 0x0625DF73;
        unsigned int y = 0xD1B84B45;
        unsigned int z = 0x152AD8D0;
        float3 f = float3(0, 0, 0);
        for(int i=0; i<32; i++)
        {
            float3 s( x/(float)0xffffffff,
                     y/(float)0xffffffff,
                     z/(float)0xffffffff);
            f += s * cos(s.dot(r));
            x = x << 1 | x >> 31;
            y = y << 1 | y >> 31;
            z = z << 1 | z >> 31;
        }
        return f * (1.0 / 64.0);
    }
};

class Diffuse : public Material
{
    float3 kd;
public:
    Diffuse(float3 kd):kd(kd) {};
    float3 shade(
                 float3 normal,
                 float3 viewDir,
                 float3 lightDir,
                 float3 lightPowerDensity,
                 float3 position)
    {
        float cosTheta = normal.dot(lightDir);
        if(cosTheta < 0) return float3(0,0,0);
        return kd * lightPowerDensity * cosTheta;
    }
};

class PhongBlinn : public Material {
    float3 ks;
    float3 kd;
    float shininess;
public:
    PhongBlinn(float3 color) {
        kd = color;
        ks = float3(0.5,0.5,0.5);
        shininess = 20.0;
    }
    float3 shade( float3 normal, float3 viewDir,
                 float3 lightDir, float3 lightPowerDensity, float3 position)
    {
        float3 halfway =(viewDir + lightDir).normalize();
        float cosDelta = normal.dot(halfway);
        if(cosDelta < 0) return float3(0,0,0);
        
        float cosTheta = normal.dot(lightDir);
        if(cosTheta < 0) return float3(0,0,0);
        
        return (lightPowerDensity * ks * pow(cosDelta, shininess)) + (kd * lightPowerDensity * cosTheta);
    }
};

class BumpyPhongBlinn : public Material {
    float3 ks;
    float3 kd;
    float shininess;
    float scale;
    float freq;
public:
    BumpyPhongBlinn(float3 color) {
        kd = color;
        ks = float3(0.5,0.5,0.5);
        shininess = 20.0;
        scale = 150;
        freq = 1.2;
    }
    float3 shade( float3 norm, float3 viewDir,
                 float3 lightDir, float3 lightPowerDensity, float3 position)
    {
        float3 normal = norm + snoiseGrad(position * scale)*freq;
        float3 halfway =(viewDir + lightDir).normalize();
        float cosDelta = normal.dot(halfway);
        if(cosDelta < 0) return float3(0,0,0);
        
        float cosTheta = normal.dot(lightDir);
        if(cosTheta < 0) return float3(0,0,0);
        
        return (lightPowerDensity * ks * pow(cosDelta, shininess)) + (kd * lightPowerDensity * cosTheta);
    }
};

float3 goldRI = float3(0.21,0.485,1.29);
float3 goldEC = float3(3.13,2.23,1.76);
float3 silverRI = float3(0.15,0.14,0.13);
float3 silverEC = float3(3.7,3.11,2.47);

class Metal : public Material {
    float3 r0;
    float3 ks;
    float shininess;
public:
    Metal(float3  refractiveIndex, float3  extinctionCoefficient){
        float3 rim1 = refractiveIndex - float3(1,1,1);
        float3 rip1 = refractiveIndex + float3(1,1,1);
        float3 k2 = extinctionCoefficient * extinctionCoefficient;
        r0 = (rim1*rim1 + k2) / (rip1*rip1 + k2);
        ks = extinctionCoefficient;
        shininess = 20;
    }
    
    float3 shade( float3 normal, float3 viewDir,
                 float3 lightDir, float3 lightPowerDensity, float3 position)
    {
        float3 halfway = (viewDir + lightDir).normalize();
        float cosDelta = normal.dot(halfway);
        if(cosDelta < 0) return float3(0,0,0);
        return (lightPowerDensity * ks * pow(cosDelta, shininess))*0.5;
    }
    
    struct Event{
        float3 reflectionDir;
        float3 reflectance;
    };
    Event evaluateEvent(float3 inDir, float3 normal) {
        Event e;
        float cosa = -normal.dot(inDir);
        float3 perp = -normal * cosa;
        float3 parallel = inDir - perp;
        e.reflectionDir = parallel - perp;
        e.reflectance = r0 + (float3(1,1,1)-r0) * pow(1 - cosa, 5);
        return e;
    }
    
};

class BumpyMetal : public Material {
    float3 r0;
    float3 ks;
    float shininess;
    float scale;
    float freq;
    
public:
    BumpyMetal(float3  refractiveIndex, float3  extinctionCoefficient){
        float3 rim1 = refractiveIndex - float3(1,1,1);
        float3 rip1 = refractiveIndex + float3(1,1,1);
        float3 k2 = extinctionCoefficient * extinctionCoefficient;
        r0 = (rim1*rim1 + k2) / (rip1*rip1 + k2);
        ks = extinctionCoefficient;
        shininess = 20;
        scale = 150;
        freq = 1.2;
    }
    
    float3 shade( float3 norm, float3 viewDir,
                 float3 lightDir, float3 lightPowerDensity, float3 position)
    {
        float3 normal = norm + snoiseGrad(position * scale)*freq;
        
        float3 halfway = (viewDir + lightDir).normalize();
        float cosTheta = normal.dot(lightDir);
        if(cosTheta < 0) return float3(0,0,0);
        float cosDelta = normal.dot(halfway);
        if(cosDelta < 0) return float3(0,0,0);
        return (lightPowerDensity * ks * pow(cosDelta, shininess))*0.5;
        
        //return float3(0,0,0);
    }
    
    struct Event{
        float3 reflectionDir;
        float3 reflectance;
    };
    Event evaluateEvent(float3 inDir, float3 norm, float3 position) {
        float3 normal = norm + snoiseGrad(position * scale)*freq;
        
        Event e;
        float cosa = -normal.dot(inDir);
        float3 perp = -normal * cosa;
        float3 parallel = inDir - perp;
        e.reflectionDir = parallel - perp;
        e.reflectance = r0 + (float3(1,1,1)-r0) * pow(1 - cosa, 5);
        return e;
    }
    
};

class Dielectric : public Material {
    float  refractiveIndex;
    float  r0;
public:
    Dielectric(float refractiveIndex): refractiveIndex(refractiveIndex) {
        r0 = (refractiveIndex - 1)*(refractiveIndex - 1)
        / (refractiveIndex + 1)*(refractiveIndex + 1);  }
    
    struct Event{
        float3 reflectionDir;
        float3 refractionDir;
        float reflectance;
        float transmittance;
    };
    
    Event evaluateEvent(float3 inDir, float3 normal) {
        Event e;
        float cosa = -normal.dot(inDir);
        float3 perp = -normal * cosa;
        float3 parallel = inDir - perp;
        e.reflectionDir = parallel - perp;
        
        float ri = refractiveIndex;
        if (cosa < 0) { cosa = -cosa; normal = -normal; ri = 1/ri; }
        float disc = 1 - (1 - cosa * cosa) / ri / ri;
        if(disc < 0) e.reflectance = 1;
        else {
            float cosb = (disc < 0)?0:sqrt(disc);
            e.refractionDir = parallel / ri - normal * cosb;
            e.reflectance = r0 + (1 - r0) * pow(1 - cosa, 5);
        }
        e.transmittance = 1 - e.reflectance;
        return e;
    }
    
};


// LIGHT SOURCES

class LightSource
{
public:
    virtual float3 getPowerDensityAt ( float3 x )=0;
    virtual float3 getLightDirAt     ( float3 x )=0;
    virtual float  getDistanceFrom   ( float3 x )=0;
};

class DirLightSource : public LightSource {
    float3 dir;
    float3 powerDensity;
public:
    DirLightSource(float3 dir, float3 powerDensity):dir(dir), powerDensity(powerDensity)
    {
    }
    float3 getPowerDensityAt( float3 x) {
        return powerDensity;
    }
    float3 getLightDirAt     ( float3 x ) {
        return dir;
    }
    float  getDistanceFrom   ( float3 x ) {
        return FLT_MAX;
    }
};

class PointLightSource : public LightSource {
    float3 position;
    float3 power;
public:
    PointLightSource(float3 position, float3 power):position(position), power(power)
    {
    }
    float3 getLightDirAt( float3 x ) {
        float3 diff = position - x;
        return diff.normalize();
    }
    float  getDistanceFrom( float3 x ) {
        float3 diff = position - x;
        return diff.norm();
    }
    float3 getPowerDensityAt( float3 x) {
        float term = 1.0/(4.0 * M_PI * getDistanceFrom(x)*getDistanceFrom(x) );
        return power*term;
    }
    
};

// MISC. SCENE ELEMENTS

// Skeletal camera class.
class Camera
{
    float3 eye;		//< world space camera position
    float3 lookAt;	//< center of window in world space
    float3 right;	//< vector from window center to window right-mid (in world space)
    float3 up;		//< vector from window center to window top-mid (in world space)
    
public:
    Camera()
    {
        //Front View
        
        eye = float3(0, 0, 3);
        lookAt = float3(0, 0, 2);
        right = float3(1, 0, 0);
        up = float3(0, 1, 0);
        
        
        
        //Birdseye View
        /*
         eye = float3(0,3,0);
         lookAt = float3(0,2,0);
         right = float3(1,0,0);
         up = float3(0,0,1);
         */
    }
    float3 getEye()
    {
        return eye;
    }
    // compute ray through pixel at normalized device coordinates
    float3 rayDirFromNdc(const float2 ndc) {
        return (lookAt - eye
                + right * ndc.x
                + up    * ndc.y
                ).normalize();
    }
};

// Ray structure.
class Ray
{
public:
    float3 origin;
    float3 dir;
    Ray(float3 o, float3 d)
    {
        origin = o;
        dir = d;
    }
};

// Hit record structure. Contains all data that describes a ray-object intersection point.
class Hit
{
public:
    Hit()
    {
        t = -1;
    }
    float t;				//< Ray paramter at intersection. Negative means no valid intersection.
    float3 position;		//< Intersection coordinates.
    float3 normal;			//< Surface normal at intersection.
    Material* material;		//< Material of intersected surface.
};

// Object abstract base class.
class Intersectable
{
protected:
    Material* material;
public:
    Intersectable(Material* material):material(material) {}
    virtual Hit intersect(const Ray& ray)=0;
};

// Simple helper class to solve quadratic equations with the Quadratic Formula [-b +- sqrt(b^2-4ac)] / 2a, and store the results.
class QuadraticRoots
{
public:
    float t1;
    float t2;
    // Solves the quadratic a*a*t + b*t + c = 0 using the Quadratic Formula [-b +- sqrt(b^2-4ac)] / 2a, and set members t1 and t2 to store the roots.
    QuadraticRoots(float a, float b, float c)
    {
        float discr = b * b - 4.0 * a * c;
        if ( discr < 0 ) // no roots
        {
            t1 = -1;
            t2 = -1;
            return;
        }
        float sqrt_discr = sqrt( discr );
        t1 = (-b + sqrt_discr)/2.0/a;
        t2 = (-b - sqrt_discr)/2.0/a;
    }
    // Returns the lesser of the positive solutions, or a negative value if there was no positive solution.
    float getLesserPositive() {
        
        return ((0 < t1 && t1 < t2) || t2 < 0)?t1:t2;
        
    }
};


// SHAPES

class Sphere : public Intersectable
{
    float3 center;
    float radius;
public:
    Sphere(const float3& center, float radius, Material* material):
    Intersectable(material),
    center(center),
    radius(radius)
    {
    }
    QuadraticRoots solveQuadratic(const Ray& ray)
    {
        float3 diff = ray.origin - center;
        float a = ray.dir.dot(ray.dir);
        float b = diff.dot(ray.dir) * 2.0;
        float c = diff.dot(diff) - radius * radius;
        return QuadraticRoots(a, b, c);
        
    }
    float3 getNormalAt(float3 r)
    {
        return (r - center).normalize();
    }
    Hit intersect(const Ray& ray)
    {
        // This is a generic intersect that works for any shape with a quadratic equation. solveQuadratic should solve the proper equation (+ ray equation) for the shape, and getNormalAt should return the proper normal
        float t = solveQuadratic(ray).getLesserPositive();
        
        Hit hit;
        hit.t = t;
        hit.material = material;
        hit.position = ray.origin + ray.dir * t;
        hit.normal = getNormalAt(hit.position);
        
        return hit;
    }
};

class Plane : public Intersectable
{
    float3 n; //plane normal
    float3 r0; //point on plane
public:
    Plane(float3 n, float3 r0, Material* material):
    Intersectable(material),
    n(n),
    r0(r0)
    {
    }
    float3 normal() {
        return n;
    }
    float3 point() {
        return r0;
    }
    Hit intersect(const Ray& ray) {
        float3 num = (r0 - ray.origin);
        float t = num.dot(n) / ray.dir.dot(n);
        
        Hit hit;
        hit.t = t;
        hit.material = material;
        hit.position = ray.origin + ray.dir * t;
        hit.normal = n;
        
        return hit;
    }
    
    
};

class Quadric : public Intersectable
{
    float4x4 coeffs = float4x4(1, 0, 0, 0,
                               0, 1, 0, 0,
                               0, 0, 1, 0,
                               0, 0, 0, -1);
public:
    Quadric(Material* material):
    Intersectable(material)
    {
    }
    Quadric* sphere() {
        coeffs = float4x4(1, 0, 0, 0,
                          0, 1, 0, 0,
                          0, 0, 1, 0,
                          0, 0, 0, -1);
        return this;
    }
    Quadric* cyl() {
        coeffs = float4x4(1, 0, 0, 0,
                          0, 0, 0, 0,
                          0, 0, 1, 0,
                          0, 0, 0, -1);
        return this;
    }
    Quadric* cone() {
        coeffs = float4x4(1, 0, 0, 0,
                          0, -1, 0, 0,
                          0, 0, 1, 0,
                          0, 0, 0, 0);
        return this;
    }
    Quadric* parab() {
        coeffs = float4x4(1, 0, 0, 0,
                          0, 0, 0, -1,
                          0, 0, 1, 0,
                          0, 0, 0, 0);
        return this;
        
    }
    Quadric* hyperb() {
        coeffs = float4x4(1, 0, 0, 0,
                          0, -1, 0, 0,
                          0, 0, 1, 0,
                          0, 0, 0, -1);
        return this;
    }
    Quadric* hparab() {
        coeffs = float4x4(1, 0, 0, 0,
                          0, 0, 0, -1,
                          0, 0, -1, 0,
                          0, 0, 0, 0);
        return this;
    }
    Quadric* hcyl() {
        coeffs = float4x4(-1, 0, 0, 0,
                          0, 0, 0, 0,
                          0, 0, 1, 0,
                          0, 0, 0, -1);
        return this;
    }
    Quadric* parallelPlanes(float height) {
        coeffs = float4x4::identity();
        coeffs._00 = 0;
        coeffs._11 = 1;
        coeffs._22 = 0;
        coeffs._33 = -(height/2)*(height/2);
        return this;
    }
    Quadric* transform(float4x4 tMatrix) {
        float4x4 tInverse = tMatrix.invert();
        float4x4 tTrans = tInverse.transpose();
        coeffs = (tInverse * coeffs) * tTrans;
        return this;
    }
    
    QuadraticRoots solveQuadratic(const Ray& ray)
    {
        float4 d = float4(ray.dir);
        d.w = 0;
        float4 e = float4(ray.origin);
        
        float a = d.dot(coeffs * d);
        float b = d.dot(coeffs * e) + e.dot(coeffs * d);
        float c = e.dot(coeffs * e);
        return QuadraticRoots(a, b, c);
        
    }
    bool contains(float3 r)
    {
        float4 rhomo(r);
        // evaluate implicit eq
        float inside = rhomo.dot(coeffs * rhomo);
        if (inside < 0) {
            return true;
        }
        else {
            return false;
        }
    }
    float3 getNormalAt(float3 r)
    {
        float4 x = float4(r);
        float4 result = operator*(x, coeffs) + coeffs.operator*(x);
        return float3(result.x,result.y,result.z).normalize();
    }
    Hit intersect(const Ray& ray)
    {
        // This is a generic intersect that works for any shape with a quadratic equation. solveQuadratic should solve the proper equation (+ ray equation) for the shape, and getNormalAt should return the proper normal
        float t = solveQuadratic(ray).getLesserPositive();
        
        Hit hit;
        hit.t = t;
        hit.material = material;
        hit.position = ray.origin + ray.dir * t;
        hit.normal = getNormalAt(hit.position);
        
        return hit;
    }
};

class ClippedQuadric: public Intersectable
{
    Quadric* shape;
    Quadric* clipper;
    
public:
    ClippedQuadric(Material* material):
    Intersectable(material)
    {
    }
    ClippedQuadric* transform(float4x4 tMatrix) {
        shape->transform(tMatrix);
        clipper->transform(tMatrix);
        return this;
    }
    ClippedQuadric* shapetransform(float4x4 tMatrix) {
        shape->transform(tMatrix);
        return this;
    }
    ClippedQuadric* cliptransform(float4x4 tMatrix) {
        clipper->transform(tMatrix);
        return this;
    }
    ClippedQuadric* sphere(float height) {
        shape = (new Quadric(material))->sphere();
        clipper = (new Quadric(material))->parallelPlanes(height);
        return this;
    }
    ClippedQuadric* cyl(float height) {
        shape = (new Quadric(material))->cyl();
        clipper = (new Quadric(material))->parallelPlanes(height);
        return this;
    }
    ClippedQuadric* cone(float height) {
        shape = (new Quadric(material))->cone();
        clipper = (new Quadric(material))->parallelPlanes(height);
        return this;
    }
    ClippedQuadric* hyperb(float height) {
        shape = (new Quadric(material))->hyperb();
        clipper = (new Quadric(material))->parallelPlanes(height);
        return this;
    }
    ClippedQuadric* parab(float height) {
        shape = (new Quadric(material))->parab();
        clipper = (new Quadric(material))->parallelPlanes(height);
        return this;
    }
    Hit intersect(const Ray& ray)
    {
        QuadraticRoots roots = shape->solveQuadratic(ray);
        
        if(roots.t1 >= 0) {
            float3 intersect1 = ray.origin + ray.dir * roots.t1;
            if (!(clipper->contains(intersect1))) {
                roots.t1 = -1;
                
            }
            
        }
        
        if(roots.t2 >= 0) {
            float3 intersect2 = ray.origin + ray.dir * roots.t2;
            if (!(clipper->contains(intersect2))) {
                roots.t2 = -1;
                
            }
        }
        
        float t = roots.getLesserPositive();
        
        Hit hit;
        hit.t = t;
        hit.material = material;
        hit.position = ray.origin + ray.dir * t;
        hit.normal = shape->getNormalAt(hit.position);
        
        return hit;
    }
    
};

class Tile: public Intersectable
{
    Plane* plane;
    float x1; // lower x bound (inclusive)
    float x2; // upper x bound (inclusive)
    float z1; // lower z bound (inclusive)
    float z2; // upper z bound (inclusive)
public:
    Tile(Material* material):
    Intersectable(material)
    {
    }
    
    //untransformable
    
    Tile* square(float xi, float xf, float zi, float zf) {
        plane = (new Plane(float3(0.0,3.0,0.0), float3(0.0,-1.0,0.0), material));
        x1 = xi;
        x2 = xf;
        z1 = zi;
        z2 = zf;
        return this;
    }
    Hit intersect(const Ray& ray) {
        float3 n = plane->normal();
        float3 r0 = plane->point();
        float3 num = (r0 - ray.origin);
        float t = num.dot(n) / ray.dir.dot(n);
        float3 position = ray.origin + ray.dir * t;
        
        if(t >= 0) {
            if (position.x < x1 || position.x > x2 || position.z < z1 || position.z > z2) {
                t = -1;
            }
        }
        
        Hit hit;
        hit.t = t;
        hit.material = material;
        hit.position = position;
        hit.normal = n;
        
        return hit;
    }
    
};

// SCENE GLOBALS

float3 background = float3(0.8,0.8,0.8);
float ep = 0.05;
int maxDepth = 2;
float3 attenuation = float3(0.001,0.001,0.001);
float3 source = float3(0,0,0);

// SCENE

class Scene
{
    Camera camera;
    std::vector<Intersectable*> objects;
    std::vector<Material*> materials;
    std::vector<LightSource*> lightsources;
public:
    Scene()
    {
        lightsources.push_back( new DirLightSource(float3(0,0.3,1),float3(1,1,1)));
        lightsources.push_back( new PointLightSource(float3(1,1,1),float3(10,10,10)));
        lightsources.push_back( new PointLightSource(float3(-1,1,1),float3(10,10,10)));
        
        materials.push_back( new BumpyPhongBlinn(float3(9.4/100, 51.8/100, 78.6/100)));// Blue Plastic at 0;
        materials.push_back( new Diffuse(float3(56.9/100, 6.3/100, 26.7/100))); //Burgundy at 1
        materials.push_back( new Diffuse(float3(12.9/100, 36.5/100, 54.5/100))); //Blue at 2
        materials.push_back( new Metal(goldRI,goldEC)); //Gold at 3
        materials.push_back( new Metal(silverRI,silverEC)); //Silver at 4
        materials.push_back( new PhongBlinn(float3(70.4/100, 14.9/100, 44.9/100)));// Red Plastic at 5
        materials.push_back( new Dielectric(1.46)); // Glass at 6
        materials.push_back( new BumpyMetal(goldRI,goldEC)); //Bumpy Gold at 7
        materials.push_back( new BumpyMetal(silverRI,silverEC)); //Bumpy Silver at 7
        
        //Chessboard
        int tilecolor = 1;
        for (float i = -2; i<2; i = i + 0.5){
            for (float j = -3; j<1; j = j + 0.5) {
                objects.push_back(( new Tile(materials.at(tilecolor)))->square(i,i+0.5,j,j+0.5));
                if (tilecolor == 1) {
                    tilecolor = 2;
                }
                else {
                    tilecolor = 1;
                }
            }
            if (tilecolor == 1) {
                tilecolor = 2;
            }
            else {
                tilecolor = 1;
            }
        }
        
        //Pawn at square 2
        objects.push_back(( new ClippedQuadric(materials.at(5)))->cone(1)->shapetransform(float4x4::translation(float3(-10,0,6))*
                                                                                          float4x4::scaling(float3(0.125,1,0.125)))->cliptransform(float4x4::translation(float3(0,-3.5,0))*float4x4::scaling(float3(1,0.25,1))));
        objects.push_back(( new ClippedQuadric(materials.at(5)))->cone(1)->shapetransform(float4x4::translation(float3(-10,0,6))*
                                                                                          float4x4::scaling(float3(0.125,1,0.125)))->cliptransform(float4x4::translation(float3(0,-1.1,0))*float4x4::scaling(float3(1,0.5,1))));
        objects.push_back(( new Quadric(materials.at(5)))->sphere()->transform(float4x4::translation(float3(-10,-2,6))*float4x4::scaling(float3(0.125,0.125,0.125))));
       
        
        //Bishop at square 3
        objects.push_back(( new ClippedQuadric(materials.at(7)))->hyperb(1)->shapetransform(float4x4::translation(float3(-15,-0.5,16))*
                                                                                            float4x4::scaling(float3(0.05,0.5,0.05)))->cliptransform(float4x4::translation(float3(0,-0.5,0))));
        objects.push_back(( new Quadric(materials.at(8)))->sphere()->transform(float4x4::translation(float3(-6,0,6.8))*float4x4::scaling(float3(0.125,0.2,0.125))));
        objects.push_back(( new ClippedQuadric(materials.at(7)))->sphere(1)->shapetransform(float4x4::translation(float3(-3.75,-5,4))*
                                                                                            float4x4::scaling(float3(0.2,0.2,0.2)))->cliptransform(float4x4::translation(float3(0,-0.5,0))));
         
        //King at square 4
        objects.push_back(( new ClippedQuadric(materials.at(6)))->hyperb(1)->shapetransform(float4x4::translation(float3(-2.5,-0.5,8))*
                                                                                            float4x4::scaling(float3(0.1,0.5,0.1)))->cliptransform(float4x4::translation(float3(0,-0.5,0))));
        objects.push_back(( new Quadric(materials.at(6)))->sphere()->transform(float4x4::translation(float3(-1,0,3))*float4x4::scaling(float3(0.25,0.25,0.25))));
        
        objects.push_back(( new ClippedQuadric(materials.at(3)))->parab(1)->shapetransform(float4x4::translation(float3(-1,0.20,3))*float4x4::scaling(float3(0.25,1,0.25))));
        
        objects.push_back(( new ClippedQuadric(materials.at(4)))->parab(1)->shapetransform(float4x4::translation(float3(-1,0.30,3))*float4x4::scaling(float3(0.25,1,0.25)))->cliptransform(float4x4::translation(float3(0,0.1,0))));
        
        objects.push_back(( new ClippedQuadric(materials.at(3)))->parab(1)->shapetransform(float4x4::translation(float3(-1,0.40,3))*float4x4::scaling(float3(0.25,1,0.25)))->cliptransform(float4x4::translation(float3(0,0.2,0))));
        
        
        //Queen at square 5
        objects.push_back(( new ClippedQuadric(materials.at(3)))->hyperb(1)->shapetransform(float4x4::translation(float3(2.5,-0.5,8))*
                                                                                            float4x4::scaling(float3(0.1,0.5,0.1)))->cliptransform(float4x4::translation(float3(0,-0.5,0))));
        objects.push_back(( new Quadric(materials.at(3)))->sphere()->transform(float4x4::translation(float3(1,0,3))*float4x4::scaling(float3(0.25,0.25,0.25))));
        for ( float t = 0.0; t < 6.29; t = t + 6.29/10) {
            objects.push_back( new Sphere(float3(cos(t)/4+0.26,0.20,sin(t)/4+0.70), 0.05, materials.at(4)) );
        }
        
                //Knight at square 6
        objects.push_back(( new ClippedQuadric(materials.at(0)))->cyl(1)->shapetransform(float4x4::translation(float3(7.5,-0.5,8))*
                                                                                         float4x4::scaling(float3(0.1,0.5,0.1)))->cliptransform(float4x4::translation(float3(0,-0.5,0))));
        objects.push_back(( new ClippedQuadric(materials.at(0)))->sphere(1)->shapetransform(float4x4::translation(float3(3.75,-5,4))*
                                                                                            float4x4::scaling(float3(0.2,0.2,0.2)))->cliptransform(float4x4::translation(float3(0,-0.5,0))));
        objects.push_back(( new Quadric(materials.at(0)))->sphere()->transform(float4x4::translation(float3(3.75,-0.5,4))*float4x4::scaling(float3(0.2,0.2,0.2))));
        objects.push_back(( new Quadric(materials.at(0)))->sphere()->transform(float4x4::translation(float3(8,-0.8,6))*float4x4::scaling(float3(0.125,0.125,0.125))));
        
    }
    
    ~Scene()
    {
        // UNCOMMENT THESE WHEN APPROPRIATE
        //for (std::vector<Material*>::iterator iMaterial = materials.begin(); iMaterial != materials.end(); ++iMaterial)
        //	delete *iMaterial;
        //for (std::vector<Intersectable*>::iterator iObject = objects.begin(); iObject != objects.end(); ++iObject)
        //	delete *iObject;
    }
    
public:
    Camera& getCamera()
    {
        return camera;
    }
    
    Hit firstIntersect(Ray ray) {
        
        Hit bestHit;
        float t = 1000000;
        for(int i = 0; i<objects.size(); i++) {
            Hit h = objects.at(i)->intersect(ray);
            if (h.t < t && h.t>0) {
                bestHit = h;
                t = h.t;
            }
        }
        return bestHit;
    }
    
    
    float3 trace(const Ray& ray, int depth, float3 s, float3 q) {
        if(depth > maxDepth) {return background;}
        Hit hit = firstIntersect(ray);
        // hit provides x, n, material
        
        // if nothing hit, return sky color
        if(hit.t < 0) return background;
        if(-ray.dir.dot(hit.normal) < 0) hit.normal = - hit.normal;
        
        float3 outRadiance = float3(0, 0, 0);
        
        Metal* metal = dynamic_cast<Metal*>(hit.material);
        if(metal != NULL){
            Metal::Event e = metal->evaluateEvent(ray.dir,hit.normal);
            outRadiance += trace(Ray(hit.position + (hit.normal * ep), e.reflectionDir), depth + 1, s, q) * e.reflectance;
        }
        
        BumpyMetal* bmetal = dynamic_cast<BumpyMetal*>(hit.material);
        if(bmetal != NULL){
            BumpyMetal::Event e = bmetal->evaluateEvent(ray.dir,hit.normal, hit.position);
            outRadiance += trace(Ray(hit.position + (hit.normal * ep), e.reflectionDir), depth + 1, s, q) * e.reflectance;
        }
        
        Dielectric* dielectric=dynamic_cast<Dielectric*>(hit.material);
        if(dielectric != NULL) {
            Dielectric::Event e = dielectric->evaluateEvent(ray.dir,hit.normal);
            outRadiance += trace(Ray(hit.position + (hit.normal * ep), e.reflectionDir), depth + 1, s, q) * e.reflectance;
            if(e.transmittance > 0)
                outRadiance += trace(Ray(hit.position - (hit.normal * ep), e.refractionDir), depth + 1, s, q) * e.transmittance;
            
            // all operations elementwise (RGB)
            outRadiance.x *= exp(- s.x * hit.t);
            outRadiance.y *= exp(- s.y * hit.t);
            outRadiance.z *= exp(- s.z * hit.t);
            outRadiance.x += q.x * (1 - exp(- s.x * hit.t)) / s.x;
            outRadiance.y += q.y * (1 - exp(- s.y * hit.t)) / s.y;
            outRadiance.z += q.z * (1 - exp(- s.z * hit.t)) / s.z;
            
        }
        
        for(unsigned int i = 0; i < lightsources.size(); i++){
            // light source provides li , Mi , |x - yi|
            
            Ray shadowRay(hit.position + (hit.normal * ep), lightsources.at(i)->getLightDirAt(hit.position));
            Hit shadowHit = firstIntersect(shadowRay);
            
            if(shadowHit.t > 0 && shadowHit.t < lightsources.at(i)->getDistanceFrom(hit.position)) continue;
            
            outRadiance += hit.material->shade(hit.normal,-ray.dir,lightsources.at(i)->getLightDirAt(hit.position), lightsources.at(i)->getPowerDensityAt(hit.position), hit.position);
            
        }
        
        return outRadiance;
        
    }
    
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
// global application data

// screen resolution
const int screenWidth = 600;
const int screenHeight = 600;
// image to be computed by ray tracing
float3 image[screenWidth*screenHeight];

Scene scene;

bool computeImage()
{
    static unsigned int iPart = 0;
    
    if(iPart >= 64)
        return false;
    for(int j = iPart; j < screenHeight; j+=64)
    {
        for(int i = 0; i < screenWidth; i++)
        {
            float3 pixelColor = float3(0, 0, 0);
            float2 ndcPixelCentre( (2.0 * i - screenWidth) / screenWidth, (2.0 * j - screenHeight) / screenHeight );
            
            Camera& camera = scene.getCamera();
            Ray ray = Ray(camera.getEye(), camera.rayDirFromNdc(ndcPixelCentre));
            
            image[j*screenWidth + i] = scene.trace(ray, 0, attenuation, source);
        }
    }
    iPart++;
    return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// OpenGL starts here. In the ray tracing example, OpenGL just outputs the image computed to the array.

// display callback invoked when window needs to be redrawn
void onDisplay( ) {
    glClearColor(0.1f, 0.2f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // clear screen
    
    if(computeImage())
        glutPostRedisplay();
    glDrawPixels(screenWidth, screenHeight, GL_RGB, GL_FLOAT, image);
    
    glutSwapBuffers(); // drawing finished
}

int main(int argc, char **argv) {
    glutInit(&argc, argv);						// initialize GLUT
    glutInitWindowSize(screenWidth, screenHeight);				// startup window size 
    glutInitWindowPosition(100, 100);           // where to put window on screen
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);    // 8 bit R,G,B,A + double buffer + depth buffer
    
    glutCreateWindow("Ray caster");				// application window is created and displayed
    
    glViewport(0, 0, screenWidth, screenHeight);
    
    glutDisplayFunc(onDisplay);					// register callback
    
    glutMainLoop();								// launch event handling loop
    
    return 0;
}