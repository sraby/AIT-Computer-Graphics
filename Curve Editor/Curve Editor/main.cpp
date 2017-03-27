//
//  main.cpp
//  Curve Editor
//
//  Created by Sam Raby on 10/2/15.
//
//





class Curve
{
public:
    // return curve point at parameter t in [0,1]
    virtual float2 getPoint(float t)=0;
    
    void draw()
    {
        glBegin(GL_LINE_STRIP);
        for (float t = 0; t <= 1; t = t + 0.01) {
            glVertex2d(getPoint(t).x,getPoint(t).y);
        }
        glEnd();
    }
};

void onDisplay() {
    
    glClearColor(0.9, 0.9, 0.9, 0.9);           // clear color
    glClear(GL_COLOR_BUFFER_BIT); // clear screen
    
    glutSwapBuffers();
    
}

void OnMouseClick(int button, int state, int x, int y)
{
    if (button == GLUT_MIDDLE_BUTTON &amp;&amp; state == GLUT_DOWN)
    {
        //store the x,y value where the click happened
        puts("Middle button clicked");
    }
}

int main(int argc, char *argv[]) {
    glutInit(&argc, argv);
    glutInitWindowSize(640, 480);
    glutInitWindowPosition(100, 100);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutCreateWindow("l01_Triangle");
    glutDisplayFunc(onDisplay);
    glutMouseFunc(OnMouseClick);
    glutIdleFunc(onIdle);
    glutMainLoop();
    
    // Event loop
    return 0;
}

