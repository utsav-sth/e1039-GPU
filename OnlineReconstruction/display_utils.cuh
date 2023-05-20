#include <helper_gl.h>
#include <GL/freeglut.h>
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>
#include <gl_menu.h>

const unsigned int window_width  = 600;
const unsigned int window_height = 500;

const unsigned int Nvars = 11;
const unsigned int Nbins = 128;

const float x0max = 150.0;
const float invpmin = 0.01;
const float invpmax = 0.20;

int win_[6];

#define LIMIT_X 0.2      // limit
#define CENTER_X (LIMIT_X-LIMIT_X)      // limit
#define LIMIT_Y 3500.0      // limit
#define CENTER_Y (LIMIT_Y-LIMIT_Y-1)      // limit

GLuint vbo;
//struct cudaGraphicsResource *gpu_tracker_display;

void setup();             // initialization of the program
//void display();           // drawing method
void createCoordinate();  // certasian coordinate
void draw();              // draw the object
void createBox(float, float, float, float);
void createHistogram(int, float*);

////////////////////////////////////////////////////////////////////////////////
//! Create VBO
////////////////////////////////////////////////////////////////////////////////
void createVBO(GLuint *vbo, float* values)
{
    assert(vbo);

    // create buffer object
    glGenBuffers(1, vbo);
    glBindBuffer(GL_ARRAY_BUFFER, *vbo);
	
    // initialize buffer object
    unsigned int size = Nbins * Nvars * sizeof(float);
    glBufferData(GL_ARRAY_BUFFER, size, values, GL_DYNAMIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, 0);

    SDK_CHECK_ERROR_GL();
}

////////////////////////////////////////////////////////////////////////////////
//! Delete VBO
////////////////////////////////////////////////////////////////////////////////
void deleteVBO(GLuint *vbo)
{

    glBindBuffer(1, *vbo);
    glDeleteBuffers(1, vbo);

    *vbo = 0;
}


void createHistogram(int npts, const float xmin, const float xmax, float *data, const int var){
    float width = (xmax-xmin)/npts;
    int i;
    for(i=0; i<npts; i++){
//	printf("  i %d x %1.4f data %1.1f  ", i, xmin+width*(i+0.5f), data[var*Nbins+i]);
	createBox(xmin+width*(i+0.5f), 1.0f, width, data[var*Nbins+i]);
    }
    printf("\n");
}

void createBox(float x, float y, float width, float height){
    glBegin(GL_POLYGON);
        glVertex2f(x, y);
        glVertex2f(x, y+height);
        glVertex2f(x+width, y+height);
        glVertex2f(x+width, y);
    glEnd();
}



void setup(){
    glClearColor(1.0, 1.0, 1.0, 1.0);
    gluOrtho2D(CENTER_X, LIMIT_X, CENTER_Y, LIMIT_Y); // -x1, x2, -y1, y2
}


void draw_data(){
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
        float* data = (float*) glMapBuffer(GL_ARRAY_BUFFER, GL_READ_ONLY);

//	for(int i=0; i<Nbins; i++)cout << " i " << i << " " << data[0*Nbins+i] << endl;
        
	//createHistogram(Nbins, -x0max, x0max, data, 0);
	createHistogram(Nbins, invpmin, invpmax, data, 2);
}

void display(){
    glClear(GL_COLOR_BUFFER_BIT);
    glColor3f(1.0, 0.0, 0.0);

    draw_data();
    createCoordinate();

    glFlush();
}

void createCoordinate(){
    glBegin(GL_LINES);
        // horizontal lines
        glVertex2f(-LIMIT_X, 0.0);
        glVertex2f(LIMIT_X, 0.0);

        // vertical lines
        glVertex2f(-LIMIT_X, -LIMIT_Y);
        glVertex2f(-LIMIT_X, LIMIT_Y);
    glEnd();
}

//bool runDisplay(int argc, char **argv, float* xpts, float* values)//, char *ref_file)
bool runDisplay(int argc, char **argv, float* values)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(window_width, window_height);
	glutCreateWindow("OR GPU display");
	glutDisplayFunc(display);
//	for(int i = 0; i<6; i++){
//	win_[i] = glutCreateSubWindow(i, window_width/3.f*(i%3), window_height/2.f*(i/3), window_width/3.f, window_height/2.f);
//	printf("window %d \n",  glutGetWindow());
//	glutSetWindow(win_[i]);
//	glutSetWindowTitle(wintitle[i].c_str());
//	}

	createVBO(&vbo, values);	
	
	setup();
//  glutInitWindowSize( 500, 500 );       /* A x A pixel screen window  */
//  glutInitWindowPosition( 100, 100 );
//  glutInitDisplayMode( GLUT_RGB | GLUT_SINGLE);
//  glutCreateWindow("Menu"); /* window title                   */
  
//  myinit();  
//  glutDisplayFunc(display_menu);         /* tell OpenGL main loop what     */
                        			/* set attributes                 */
//  glutMouseFunc(mouse_func);		/* Set actions on mouse Click 	*/
//  glutPassiveMotionFunc(hover);		/* Set actions on mouse hover 	*/

	glutMainLoop();
	
	return true;
}










