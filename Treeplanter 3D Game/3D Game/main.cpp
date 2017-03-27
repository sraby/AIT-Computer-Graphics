#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
// Needed on MsWindows
#include <windows.h>
#endif // Win32 platform

// Download glut from: http://www.opengl.org/resources/libraries/glut/
#include <GLUT/glut.h>
#include "float2.h"
#include "float3.h"
#include <vector>
#include <map>
#include "Mesh.h"
#include <iostream>
#include <string.h>

extern "C" unsigned char* stbi_load(char const *filename, int *x, int *y, int *comp, int req_comp);

// SET TEXTURE/OBJECT PATH NAMES HERE

const char* tiggerobj = "/Users/Sam/AIT Computer Graphics/AitGraphics-Hello/Projects/ObjectFiles/tigger.obj"; 
const char* treeobj = "/Users/Sam/AIT Computer Graphics/AitGraphics-Hello/Projects/ObjectFiles/tree.obj";
const char* tiggertexture = "/Users/Sam/AIT Computer Graphics/AitGraphics-Hello/Projects/3D Game/3D Game/mesh/tigger.png";
const char* flowerstexture = "/Users/Sam/AIT Computer Graphics/AitGraphics-Hello/Projects/ObjectFiles/flowers.jpg";
const char* tree1texture = "/Users/Sam/AIT Computer Graphics/AitGraphics-Hello/Projects/ObjectFiles/tree1.png";
const char* watertexture = "/Users/Sam/AIT Computer Graphics/AitGraphics-Hello/Projects/ObjectFiles/water.jpg";

// HELPER FNS

float length(float3 start, float3 end) {
    
    float dx = end.x - start.x;
    float dy = end.y - start.y;
    float dz = end.z - start.z;
    float x = dx * dx;
    float y = dy * dy;
    float z = dz * dz;
    
    return sqrtf(x + y + z);
}

// LIGHTING

class LightSource
{
public:
    virtual float3 getRadianceAt  ( float3 x )=0;
    virtual float3 getLightDirAt  ( float3 x )=0;
    virtual float  getDistanceFrom( float3 x )=0;
    virtual void   apply( GLenum openglLightName )=0;
};

class DirectionalLight : public LightSource
{
    float3 dir;
    float3 radiance;
public:
    DirectionalLight(float3 dir, float3 radiance)
    :dir(dir), radiance(radiance){}
    float3 getRadianceAt  ( float3 x ){return radiance;}
    float3 getLightDirAt  ( float3 x ){return dir;}
    float  getDistanceFrom( float3 x ){return 900000000;}
    void   apply( GLenum openglLightName )
    {
        float aglPos[] = {dir.x, dir.y, dir.z, 0.0f};
        glLightfv(openglLightName, GL_POSITION, aglPos);
        float aglZero[] = {0.0f, 0.0f, 0.0f, 0.0f};
        glLightfv(openglLightName, GL_AMBIENT, aglZero);
        float aglIntensity[] = {radiance.x, radiance.y, radiance.z, 1.0f};
        glLightfv(openglLightName, GL_DIFFUSE, aglIntensity);
        glLightfv(openglLightName, GL_SPECULAR, aglIntensity);
        glLightf(openglLightName, GL_CONSTANT_ATTENUATION, 1.0f);
        glLightf(openglLightName, GL_LINEAR_ATTENUATION, 0.0f);
        glLightf(openglLightName, GL_QUADRATIC_ATTENUATION, 0.0f);
    }
};

class PointLight : public LightSource
{
    float3 pos;
    float3 power;
public:
    PointLight(float3 pos, float3 power)
    :pos(pos), power(power){}
    float3 getRadianceAt  ( float3 x ){return power*(1/(x-pos).norm2()*4*3.14);}
    float3 getLightDirAt  ( float3 x ){return (pos-x).normalize();}
    float  getDistanceFrom( float3 x ){return (pos-x).norm();}
    void   apply( GLenum openglLightName )
    {
        float aglPos[] = {pos.x, pos.y, pos.z, 1.0f};
        glLightfv(openglLightName, GL_POSITION, aglPos);
        float aglZero[] = {0.0f, 0.0f, 0.0f, 0.0f};
        glLightfv(openglLightName, GL_AMBIENT, aglZero);
        float aglIntensity[] = {power.x, power.y, power.z, 1.0f};
        glLightfv(openglLightName, GL_DIFFUSE, aglIntensity);
        glLightfv(openglLightName, GL_SPECULAR, aglIntensity);
        glLightf(openglLightName, GL_CONSTANT_ATTENUATION, 0.0f);
        glLightf(openglLightName, GL_LINEAR_ATTENUATION, 0.0f);
        glLightf(openglLightName, GL_QUADRATIC_ATTENUATION, 0.25f / 3.14f);
    }
};

// MATERIALS/TEXTURES

class Material
{
public:
    float3 kd;			// diffuse reflection coefficient
    float3 ks;			// specular reflection coefficient
    float shininess;	// specular exponent
    Material()
    {
        kd = float3(0.5, 0.5, 0.5) + float3::random() * 0.5;
        ks = float3(1, 1, 1);
        shininess = 15;
    }
    void defcolor(float3 color)
    {
        kd = color;
    }
    virtual void apply()
    {
        glDisable(GL_TEXTURE_2D);
        float aglDiffuse[] = {kd.x, kd.y, kd.z, 1.0f};
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, aglDiffuse);
        float aglSpecular[] = {kd.x, kd.y, kd.z, 1.0f};
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, aglSpecular);
        if(shininess <= 128)
            glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shininess);
        else
            glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 128.0f);
    }
};

class TexturedMaterial : public Material
{
public:
    GLuint textureName;
    
    TexturedMaterial(const char* filename,
                     GLuint filtering = GL_LINEAR_MIPMAP_LINEAR
                     ){
        unsigned char* data;
        int width;
        int height;
        int nComponents = 4;
        
        data = stbi_load(filename, &width, &height, &nComponents, 0);
        
        if(data == NULL) return;
        
        // opengl texture creation comes here
        
        glGenTextures(1, &textureName);  // id generation
        glBindTexture(GL_TEXTURE_2D, textureName);      // binding
        
        /* if(nComponents == 4)
            
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, data); //uploading
        
        else if(nComponents == 3)
            
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data); //uploading
        */
       if(nComponents == 4)
            gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGBA, width, height, GL_RGBA, GL_UNSIGNED_BYTE, data);
        else if(nComponents == 3)
            gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGB, width, height, GL_RGB, GL_UNSIGNED_BYTE, data);

        delete data;
    }
    
    void apply() {
        
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        glEnable(GL_TEXTURE_2D);
        Material::apply();
        glTexParameteri(GL_TEXTURE_2D,
                        GL_TEXTURE_MIN_FILTER,GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D,
                        GL_TEXTURE_MAG_FILTER,GL_LINEAR_MIPMAP_LINEAR);
        glTexEnvi(GL_TEXTURE_ENV,
                  GL_TEXTURE_ENV_MODE, GL_REPLACE);
        glEnable(GL_TEXTURE_2D);
        glTexEnvi(GL_TEXTURE_ENV,
                  GL_TEXTURE_ENV_MODE,
                 GL_REPLACE); //GL_MODULATE
        glBindTexture(GL_TEXTURE_2D, textureName);
        glEnable(GL_LINEAR_MIPMAP_LINEAR);
        //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    }
    
};

//OBJECTS

class Object
{
protected:
    Material* material;
    float3 scaleFactor;
    float3 position;
    float3 orientationAxis;
    float orientationAngle;
public:
    int life = 8;
    int wins = 0;
    bool grow = false;
    bool erase = false;
    bool visible = true;
    Object(Material* material):material(material),orientationAngle(0.0f),scaleFactor(1.0,1.0,1.0),orientationAxis(0.0,1.0,0.0){}
    virtual ~Object(){}
    Object* translate(float3 offset){
        position += offset; return this;
    }
    Object* scale(float3 factor){
        scaleFactor *= factor; return this;
    }
    Object* rotate(float angle){
        orientationAngle += angle; return this;
    }
    virtual void draw()
    {
        material->apply();
        // apply scaling, translation and orientation
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glTranslatef(position.x, position.y, position.z);
        glRotatef(orientationAngle, orientationAxis.x, orientationAxis.y, orientationAxis.z);
        glScalef(scaleFactor.x, scaleFactor.y, scaleFactor.z);
        drawModel();
        glPopMatrix();
    }
    virtual void drawModel()=0;
    virtual void drawShadow(float3 lightDir)
    {
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        
        //glTranslatef(0, 0.01, 0);
       // glScalef(1, 0.01, 1);
        float shear[] = {
            1, 0, 0, 0,
            -lightDir.x/lightDir.y, 0, -lightDir.z/lightDir.y, 0,
            0, 0, 1, 0,
            0, 0, 0, 1 };
        glMultMatrixf(shear);
        
        glTranslatef(position.x, position.y, position.z);
        glRotatef(orientationAngle, orientationAxis.x, orientationAxis.y, orientationAxis.z);
        glScalef(scaleFactor.x, scaleFactor.y, scaleFactor.z);
        
        drawModel();
        glPopMatrix();
    }
    virtual void move(double t, double dt){}
    virtual bool control(std::vector<bool>& keysPressed, std::vector<Object*>& spawn, std::vector<Object*>& objects){return false;}
    virtual float3 getCenter() {
        return position;
    }
    virtual float getRadius() =0;
    void recolor(Material* newcolor) {
        material = newcolor;
    }
    void loselife() {
        life --;
    }
    void drawBitmapText(char *string,float x,float y,float z)
    {
        char *c;
        glRasterPos3f(x, y,z);
        
        for (c=string; *c != ' '; c++)
        {
            glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10, *c);
        }
    }
};

class Teapot : public Object
{
public:
    Teapot(Material* material):Object(material){}
    void drawModel()
    {
        glutSolidTeapot(0.5f);
    }
    float getRadius() {return 0.5;}
};


class Ground : public Object
{
public:
    Ground(Material* material):Object(material){
        position = float3(0,0,0);
    }
    void drawModel()
    {
    
        glBegin(GL_QUADS);
        
    for (double i = -200; i < 200; i = i + 100) {
            for (double j = -200; j < 200; j = j + 100)
            {
                
                glTexCoord2d(0, 0);
                glVertex3d(i,-0.01,j);
                glTexCoord2d(1, 0);
                glVertex3d(i+100,-0.01,j);
                glTexCoord2d(1, 1);
                glVertex3d(i+100,-0.01,j + 100);
                glTexCoord2d(0, 1);
                glVertex3d(i,-0.01,j + 100);
                
            }
        }
        glEnd();
    }
    void drawShadow( float3 lightDir) {
    }
    float getRadius() {return 0;}
};

class MeshInstance : public Object
{
protected:
    Mesh* mesh;
public:
    MeshInstance(Mesh* mesh, Material* material):Object(material), mesh(mesh) {}
    void drawModel()
    {
        mesh->draw();
    }
    float3 getCenter() {
        float3 sum = float3(0,0,0);
        for (unsigned int i = 0; i < mesh->positions.size(); i++) {
            sum += *mesh->positions.at(i);
        }
        float m = 1/((float)(mesh->positions.size()));
        return position + float3(sum.x * m, sum.y * m, sum.z * m);
    }
    float getRadius() {
        float3 center = getCenter();
        float radius = 0;
        for (unsigned int i = 0; i < mesh->positions.size(); i++) {
            float3 vertex = *mesh->positions.at(i);
            float len = length(center,position + vertex);
            if (len >= radius) radius = len;
        }
        
        return radius;
    }
};

// Camera Globals...
float camAngle = 0;
float3 camAvatar = float3(0,0,0);

class Bouncer : public MeshInstance
{
protected:
    float3 velocity = float3(0,10,0);
    float3 acceleration = float3(0,-10,0);
    float angularvelocity = 0;
    float angularacceleration = 0;
    float restitution = 0.95;
public:
    Bouncer(Mesh* mesh, Material* material): MeshInstance(mesh,material) {
    }
    void move(double t, double dt){
        if (wins >= 10) {
            acceleration = float3(0,2,0);
            velocity = float3(0,2,0) + acceleration*dt;
            position += velocity*dt;
            angularvelocity = 200;
            orientationAngle += angularvelocity*dt;
        }
        else if (life <= 0) {
            velocity = float3(0,0,0);
        }
        else {
        velocity += acceleration*dt;
        if(position.y < 0 && velocity.y < 0) velocity.y *= -restitution;
        if(position.x > 200 || position.x < -200 || position.z > 200 || position.z < -200)
        {
            position = float3(0,0,0);
            camAvatar = float3(0,0,0);
        }
        position += velocity*dt;
        camAvatar += velocity*dt;
        velocity.x *= pow(0.6, dt);
        velocity.z *= pow(0.6, dt);
        
        angularvelocity += angularacceleration*dt;
        orientationAngle += angularvelocity*dt;
        camAngle += (M_PI/180)*angularvelocity*dt;
        angularvelocity *= pow(0.1, dt);
        }
    }
    bool control(std::vector<bool>& keysPressed, std::vector<Object*>& spawn, std::vector<Object*>& objects){
       /* for (unsigned int iObject=0; iObject<objects.size(); iObject++)
            if (objects.at(iObject) != this && length(getCenter(), objects.at(iObject)->getCenter()) < getRadius() + objects.at(iObject)->getRadius())) {
                MeshInstance* trees = dynamic_cast<MeshInstance*>(objects.at(iObject));
                if(trees != NULL) {
                    velocity = -velocity;
                }
            } */
        
        if (keysPressed.at('a')) {
            angularacceleration += 10;
        }
        else if (keysPressed.at('d')) {
            angularacceleration += -10;
        }
        else {angularacceleration = 0;}

        if (keysPressed.at('w')) {
            acceleration = float3(-cos((M_PI/180)*orientationAngle)*10, -10, sin((M_PI/180)*orientationAngle) *10);
        }
        else if (keysPressed.at('s')) {
            acceleration = float3(cos((M_PI/180)*orientationAngle)*10, -10, -sin((M_PI/180)*orientationAngle) *10);
            
        }
        else {acceleration = float3(0,-10,0);}
        return false;
    }
};

class SpinTeapot : public Object
{
public:
    SpinTeapot(Material* material):Object(material){}
    void drawModel()
    {
        glutSolidTeapot(2.0f);
    }
    float getRadius() {return 2;}
    void move(double t, double dt) {
        orientationAngle += 50*dt;
    }
    bool control(std::vector<bool>& keysPressed, std::vector<Object*>& spawn, std::vector<Object*>& objects)
    {
        for (unsigned int iObject=0; iObject<objects.size(); iObject++) {
            if (objects.at(iObject) != this && length(getCenter(), objects.at(iObject)->getCenter()) < getRadius() + 0.85*objects.at(iObject)->getRadius()) {
                Bouncer* avatar = dynamic_cast<Bouncer*>(objects.at(iObject));
                if(avatar != NULL) {
                    visible = false;
                    position = float3(0,1000,0);
                    objects.at(iObject)->life = life + 3;
                    objects.at(iObject)->scale(float3(1,1.11,1));
                    objects.at(iObject)->scale(float3(1,1.11,1));
                    objects.at(iObject)->scale(float3(1,1.11,1));
                    
                }
            }
        }
        return false;
    }
};

class Bullet : public Teapot
{
protected:
    float3 velocity = float3(10*sin(camAngle - M_PI/2), 0, 10*cos(camAngle - M_PI/2));
    float3 acceleration = float3(0,-10,0);
    float restitution = 0.9;
    int lifetime = 250;
public:
    Bullet(Material* material):Teapot(material) {
        position = camAvatar + float3(0,3.5,0);
    }
    void move(double t, double dt){
        velocity += acceleration*dt;
        if(position.y < 0.3 && velocity.y < 0) velocity.y *= -restitution;
        position += velocity*dt;
        velocity.x *= pow(0.4, dt);
        velocity.z *= pow(0.4, dt);
        lifetime -= dt;
    }
    bool control(std::vector<bool>& keysPressed, std::vector<Object*>& spawn, std::vector<Object*>& objects) {
        if (lifetime <= 0) {
            for (int i = 0; i < objects.size(); i++) {
                if (objects.at(i) == this) visible = false;
            }
        }
        return false;
    }
};

class Dirt : public Object
{
public:
    Dirt(Material* material):Object(material) {}
    void drawModel()
    {
        glutSolidCube(3);
    }
    float getRadius() {return 3.0;}
    bool control(std::vector<bool>& keysPressed, std::vector<Object*>& spawn, std::vector<Object*>& objects) {
        for (unsigned int iObject=0; iObject<objects.size(); iObject++)
            if (objects.at(iObject) != this && length(getCenter(), objects.at(iObject)->getCenter()) < getRadius() + objects.at(iObject)->getRadius()) {
                Bullet* bullet = dynamic_cast<Bullet*>(objects.at(iObject));
                if(bullet != NULL) {
                    grow = true;
                }
            }
        return false;
    }
};

// CAMERA

class Camera
{
    
    float3 eye;
    
    float3 ahead;
    float3 lookAt;
    float3 right;
    float3 up;
    
    float fov;
    float aspect;
    
    float2 lastMousePos;
    float2 mouseDelta;
    
public:
    float3 getEye()
    {
        return eye;
    }
    Camera()
    {
        eye = float3(0, 7, -10);
        lookAt = float3(0, 0, 0);
        right = float3(1, 0, 0);
        up = float3(0, 1, 0);
        
        fov = 1.1;
        aspect  = 1;
    }
    
    void apply()
    {
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        gluPerspective(fov /M_PI*180, aspect, 0.1, 500);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        gluLookAt(eye.x, eye.y, eye.z, lookAt.x, lookAt.y, lookAt.z, 0.0, 1.0, 0.0);
    }
    
    void setAspectRatio(float ar) { aspect= ar; }
    
    void move(float dt, std::vector<bool>& keysPressed)
    {
        float3 camDirection = float3(camAvatar.x + sin(camAngle - M_PI/2), 7, camAvatar.z + cos(camAngle - M_PI/2));
        float3 camPosDir = camDirection - camAvatar;
        eye = float3(camAvatar.x - 15*camPosDir.x, 7, camAvatar.z - 15*camPosDir.z);
        lookAt = float3(camAvatar.x, 7, camAvatar.z);
        
        
        /*
        
        if(keysPressed.at('w'))
            eye += ahead * dt * 20;
        if(keysPressed.at('s'))
            eye -= ahead * dt * 20;
        if(keysPressed.at('a'))
            eye -= right * dt * 20;
        if(keysPressed.at('d'))
            eye += right * dt * 20;
        if(keysPressed.at('q'))
            eye -= float3(0,1,0) * dt * 20;
        if(keysPressed.at('e'))
            eye += float3(0,1,0) * dt * 20;
        
        float yaw = atan2f( ahead.x, ahead.z );
        float pitch = -atan2f( ahead.y, sqrtf(ahead.x * ahead.x + ahead.z * ahead.z) );
        
        yaw -= mouseDelta.x * 0.02f;
        pitch += mouseDelta.y * 0.02f;
        if(pitch > 3.14/2) pitch = 3.14/2;
        if(pitch < -3.14/2) pitch = -3.14/2;
        
        mouseDelta = float2(0, 0);
        
        ahead = float3(sin(yaw)*cos(pitch), -sin(pitch), cos(yaw)*cos(pitch) );
        right = ahead.cross(float3(0, 1, 0)).normalize();
        up = right.cross(ahead);
        
        lookAt = eye + ahead;
         
         */
        
    }
    void startDrag(int x, int y)
    {
        lastMousePos = float2(x, y);
    }
    void drag(int x, int y)
    {
        float2 mousePos(x, y);
        mouseDelta = mousePos - lastMousePos;
        lastMousePos = mousePos;
    }
    void endDrag()
    {
        mouseDelta = float2(0, 0);
    }
};

// SCENE

class Scene
{
    Camera camera;
public:
    std::vector<LightSource*> lightSources;
    std::vector<Object*> objects;
    std::vector<Material*> materials;
    std::vector<Mesh*> meshes;
    void initialize()
    {
        // BUILD YOUR SCENE HERE
        lightSources.push_back(
                               new DirectionalLight(
                                                    float3(3, 3, 3),
                                                    float3(1, 1, 1)));
        lightSources.push_back(
                               new PointLight(
                                              float3(3, 3, 3),
                                              float3(0.2, 0.1, 0.1)));
        Mesh* mesh1 = new Mesh(tiggerobj);
        meshes.push_back(mesh1);
        Mesh* mesh2 = new Mesh(treeobj);
        meshes.push_back(mesh2);
        
        TexturedMaterial* tigger = new TexturedMaterial(tiggertexture);
        materials.push_back(tigger); //tigger texture at 0
        
        TexturedMaterial* flower = new TexturedMaterial(flowerstexture);
        materials.push_back(flower); //concrete texture at 1
        
        TexturedMaterial* tree = new TexturedMaterial(tree1texture);
        materials.push_back(tree); //tree texture at 2
        
        Material* dirt = new Material();
        materials.push_back(dirt);
        dirt->kd = float3(33.3/100, 22.4/100, 0); //dirt at 3

        
        Material* redDiffuseMaterial = new Material();
        materials.push_back(redDiffuseMaterial);
        redDiffuseMaterial->kd = float3(0.8, 0.2, 0.1); //red at 4
        
        Material* yellowDiffuseMaterial = new Material();
        materials.push_back(yellowDiffuseMaterial);
        yellowDiffuseMaterial->kd = float3(51.8/100, 95.7/100, 15.7/100); //green at 5
        
        TexturedMaterial* water = new TexturedMaterial(watertexture);
        materials.push_back(water); //water texture at 6
        
        materials.push_back(new Material());
        materials.push_back(new Material());
        materials.push_back(new Material());
        materials.push_back(new Material());
        materials.push_back(new Material());
        materials.push_back(new Material());
        
        objects.push_back( (new Bouncer(meshes.at(0), tigger))->scale(float3(0.25,0.25,0.25))->translate(float3(0,0,0))); //Avatar at origin
        objects.push_back( (new Ground(flower))); //ground
        for (float i = 0; i < 6.2; i = i + M_PI/5) {
            objects.push_back( (new Dirt(materials.at(3)))->translate(float3(100*cos(i),1.5,100*sin(i))));
        } //dirt cubes
        
        for (int i = 0; i < 5; i++) {
            float3 pos = float3::random();
            pos.y = 0;
            objects.push_back( (new MeshInstance(meshes.at(1),materials.at(2)))->translate(pos * 170));
            float3 pos1 = float3::random();
            pos1.y = 0;
            pos1.x = -pos1.x;
            objects.push_back( (new MeshInstance(meshes.at(1),materials.at(2)))->translate(pos1 * 170));
            float3 pos2 = float3::random();
            pos2.y = 0;
            pos2.z = -pos2.z;
            objects.push_back( (new MeshInstance(meshes.at(1),materials.at(2)))->translate(pos2 * 170));
            float3 pos3 = float3::random();
            pos3.y = 0;
            pos3.x = -pos3.x;
            pos3.z = -pos3.z;
            objects.push_back( (new MeshInstance(meshes.at(1),materials.at(2)))->translate(pos3 * 170));
            
        }
        objects.push_back( (new SpinTeapot(materials.at(6)))->translate(float3(-25,5,5)));
        objects.push_back( (new SpinTeapot(materials.at(6)))->translate(float3(190,5,190)));
        objects.push_back( (new SpinTeapot(materials.at(6)))->translate(float3(-190,5,-190)));
    }
    
    ~Scene()
    {
        for (std::vector<LightSource*>::iterator iLightSource = lightSources.begin(); iLightSource != lightSources.end(); ++iLightSource)
            delete *iLightSource;
        for (std::vector<Material*>::iterator iMaterial = materials.begin(); iMaterial != materials.end(); ++iMaterial)
            delete *iMaterial;
        for (std::vector<Object*>::iterator iObject = objects.begin(); iObject != objects.end(); ++iObject)
            delete *iObject;
    }
    Camera& getCamera()
    {
        return camera;
    }
    
    void draw()
    {
        camera.apply();
        unsigned int iLightSource=0;

        for (; iLightSource<lightSources.size(); iLightSource++)
        {
            glEnable(GL_LIGHT0 + iLightSource);
            lightSources.at(iLightSource)->apply(GL_LIGHT0 + iLightSource);
        }
        for (; iLightSource<GL_MAX_LIGHTS; iLightSource++){
            glDisable(GL_LIGHT0 + iLightSource);
        }
        
        for (unsigned int iObject=0; iObject<objects.size(); iObject++){
            if (objects.at(iObject)->visible) objects.at(iObject)->draw();
        }
        
        float3 lightDir =
        lightSources.at(0)
        ->getLightDirAt(float3(0, 0, 0));
        
        glDisable(GL_LIGHTING);
        glDisable(GL_TEXTURE_2D);
        
        glColor3d(0.1,0.1,0.1);
        for (unsigned int iObject=0; iObject<objects.size(); iObject++){
            if (objects.at(iObject)->visible) {
            objects.at(iObject)->drawShadow(lightDir);
                }
            }
        
        glEnable(GL_LIGHTING);
    }
    void control(std::vector<bool>& keysPressed)
    {
        std::vector<Object*> spawn;
        for (unsigned int iObject=0;
             iObject<objects.size(); iObject++) {
            objects.at(iObject)->control(keysPressed, spawn, objects);
        if (objects.at(iObject)->grow) {
            objects.push_back( (new MeshInstance(meshes.at(1),materials.at(2)))->translate(objects.at(iObject)->getCenter()) );
            objects.at(iObject)->translate(float3(0,10000,0));
            objects.at(0)->wins ++;
            std::cout << objects.at(0)->wins;
        }
        }
    }
    void filter()
    {
        for (unsigned int iObject=0;
             iObject<objects.size(); iObject++) {
            if (objects.at(iObject)->grow || objects.at(iObject)->erase) {
                objects.erase(objects.begin()+iObject);
            }
        }
    }
    void move(double t, double dt)
    {
        for (unsigned int iObject=0; iObject<objects.size(); iObject++)
            objects.at(iObject)->move(t,dt);
    }
};

Scene scene;
std::vector<bool> keysPressed;
int cooldown = 0;

// INTERFACE FNS

void onDisplay( ) {
    glClearColor(0.40f, 0.69f, 0.80f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // clear screen
    
    scene.draw();
    
    glutSwapBuffers(); // drawing finished
}

void onIdle()
{
    double t = glutGet(GLUT_ELAPSED_TIME) * 0.001;        	// time elapsed since starting this program in msec
    static double lastTime = 0.0;
    double dt = t - lastTime;
    lastTime = t;
    scene.getCamera().move(dt, keysPressed);
    scene.move(t,dt);
    scene.control(keysPressed);
    scene.filter();
    if(keysPressed.at('b') && cooldown <= 0 && scene.objects.at(0)->life > 0){
        scene.objects.push_back( (new Bullet(scene.materials.at(6))));
        cooldown = 20;
        scene.objects.at(0)->scale(float3(1,0.9,1));
        scene.objects.at(0)->loselife();
    }
    if(scene.objects.at(0)->wins >= 10) {
        scene.objects.at(0)->recolor(scene.materials.at(5));
    }
    else if(scene.objects.at(0)->life <= 0) {
        scene.objects.at(0)->recolor(scene.materials.at(4));
    }
    cooldown -= 1;
    glutPostRedisplay();
}

void onKeyboard(unsigned char key, int x, int y)
{
    keysPressed.at(key) = true;
}

void onKeyboardUp(unsigned char key, int x, int y)
{
    keysPressed.at(key) = false;
}

void onMouse(int button, int state, int x, int y)
{
    if(button == GLUT_LEFT_BUTTON)
        if(state == GLUT_DOWN)
            scene.getCamera().startDrag(x, y);
        else scene.getCamera().endDrag();
}

void onMouseMotion(int x, int y)
{
    scene.getCamera().drag(x, y);
}

void onReshape(int winWidth, int winHeight)
{
    glViewport(0, 0, winWidth, winHeight);
    scene.getCamera().setAspectRatio((float)winWidth/winHeight);
}	

int main(int argc, char **argv) {
    
    glutInit(&argc, argv);						// initialize GLUT
    glutInitWindowSize(600, 600);				// startup window size 
    glutInitWindowPosition(100, 100);           // where to put window on screen
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);    // 8 bit R,G,B,A + double buffer + depth buffer
    
    glutCreateWindow("Save the Forest!");				// application window is created and displayed
    
    glViewport(0, 0, 600, 600);
    
    glutDisplayFunc(onDisplay);					// register callback
    glutIdleFunc(onIdle);						// register callback
    glutReshapeFunc(onReshape);
    glutKeyboardFunc(onKeyboard);
    glutKeyboardUpFunc(onKeyboardUp);
    glutMouseFunc(onMouse);
    glutMotionFunc(onMouseMotion);
    
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_NORMALIZE);
    
    scene.initialize();
    for(int i=0; i<256; i++)
        keysPressed.push_back(false);
    
    glutMainLoop();								// launch event handling loop
    
    return 0;
}