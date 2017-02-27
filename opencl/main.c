/*
 ===========================================================================
 Copyright (C) 2017 Jon Rood.
 
 This file is part of Enlightning source code.
 
 Enlightning source code is free software; you can redistribute it
 and/or modify it under the terms of the GNU General Public License as
 published by the Free Software Foundation; either version 3 of the License,
 or (at your option) any later version.
 
 Enlightning source code is distributed in the hope that it will be
 useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with Enlightning; if not, see <http://www.gnu.org/licenses/>.
 ===========================================================================
 */

// Notes: Using floats to allow execution on older GPUs. Calculations break down quickly.

#include <GLUT/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <sys/stat.h>
#include <OpenCL/opencl.h>

#define SQUARE(x) ((x)*(x))             //faster way to square values that using pow()
#define PI 3.14159265358979323846264338
#define KSIZE 4                         //number of equations in set of equations

//use defined simulation parameters
int X = 256;                            //size of grid in x direction
int Y = 256;                            //size of grid in y direction
int bc_left = 1;                        //left boundary condition: 0 none, 1 reflecting, 2 absorbing, 3 for periodic
int bc_right = 1;                       //right boundary condition: 0 none, 1 reflecting, 2 absorbing, 3 for periodic
int bc_bottom = 1;                      //bottom boundary condition: 0 none, 1 reflecting, 2 absorbing, 3 for periodic
int bc_top = 1;                         //top boundary condition: 0 none, 1 reflecting, 2 absorbing, 3 for periodic
int time_steps = 2000;                  //number of time steps to run simulation for
int ppw = 10;                           //number of points to use per wavelength (NEED X AND Y FOR THIS TO HAVE DIF. dy)
int frequency_max = 200;                //maximum frequency to account for
int src_x_pos = 51;                     //position of middle of point source in x
int src_y_pos = 51;                     //position of middle of point source in y
int source_mode = 0;                    //0 for pulse source, 1 for driven
int source_type = 1;                    //1 for point, 2 for sine curve, 3 for bitmap input
int paused = 0;                         //used to pause simulation
int pabsorb = 1;                        //turn on to absorb left moving wave during periodic boundary conditions
int winx = 1344;                        //window size x
int winy = 840;                         //window size y
float sine_src_amp = 8;                 //amplitude of sine wave source shape
float sine_src_frq = 0.125;             //frequency of sine wave source shape
float grid_amp = 1.0;                   //zoom z axis of grid
float grid_rot_z = 0;                   //rotate grid on z axis
float grid_rot_x = 0;                   //rotate grid on x axis
float grid_fov = 50;                    //zoom in or out on the grid
float grid_look_x = 2.0;                //move camera focus around x
float grid_look_y = 2.8;                //move camera focus around y
float cfl = 0.125;//0.49;               //<0.5 for WENO method, smaller when using molecular relaxation
float source_amp = 1e2;                 //density amplitude of source
float alph = 4;                         //wideness of source
float omega = 1000;                     //frequency of source
float alpha = 390;                      //the alpha for splitting, chosen between 300-400
float gthres = 0.00001;                 //gradient threshold for switching from DRP(1) to WENO(0)
float bmax = 0.0;                       //used for finding global pressure maximum for coloring the grid
unsigned int step = 0;                  //current step of simulation

//more advanced simulation parameters
const float c = 343;                    //speed of sound
const float gamma1 = 1.402;             //ratio of specific heats
const float p0 = 101325;                //ambient pressure in Pa
const float p_ref = 101325;             //same as p0
const float RH = 0.2;                   //relative humidity percentage
const float T0 = 293.16;                //ambient temperature in K
const float T_ref = 293.16;             //same as T0
const float rho0 = 1.21;                //ambient density
const float R = 287.06;                 //gas constant
const float R_tilde = 8314;             //another gas constant
const float M = 28.96;                  //molar mass of air
const float c_v = 720.4;                //specific heat of air at constant volume
const float c_p = 1010;                 //specific heat of air at constant pressure
const float s0 = 0;                     //ambient entropy
const float T_star_n = 3352;            //molecular constant of nitrogen
const float T_star_o = 2239;            //molecular constant of oxygen
const float mu = 1.846e-5;              //shear viscosity, 0 to turn off shear viscosity
const float muB = 0.6*1.846e-5;         //bulk viscosity, 0 to turn off bulk viscosity
const float kappa = 0.02624;            //constant of thermal conductivity, 0 to turn off thermal conductivity

//other simulation global variables
float dx;                               //width of grid element in x direction
float dy;                               //width of grid element in y direction
float dt;                               //width of time step
float x_length;                         //total length of x grid in meters
float y_length;                         //total length of y grid in meters
float lambda;                           //minimum wavelength to account for
float tau_n, tau_o;                     //apparent vibration temperature of nitrogen and oxygen
float *p;                               //pressure grid
float *T;                               //temperature grid
float *w;                               //for current time step
float *w_n;                             //for sub steps between next time step
float *K;                               //scratch grid for calculating RHS
float *b;                               //interpolated grid for opengl
float *H0;                              //source grid

//opencl objects
cl_context          context;
cl_command_queue    cmd_queue;
cl_device_id        devices;
cl_int              err;
cl_program          program;
cl_kernel           kernel;
cl_mem              w_n_mem;
cl_mem              K_mem;
cl_mem              p_mem;
cl_mem              H0_mem;

void initialConditions(void);           //sets up initial conditions and allocates some memory
void updateVariables(void);             //update all variables
void rungeKutta(void);                  //total variation diminishing runge kutta scheme
void addSource(void);                   //function that adds a smooth point source density
void bc(void);                          //apply boundary conditions
void freeMemory(void);                  //free malloc() stuff from initialConditions()
void fail(char *message);               //used to kill program if error
void mfail(char *message);              //used to kill program if malloc error
void printOutput(void);                 //prints simulation info to terminal
void readAndCheckInput(void);           //read input file and check that it makes sense
void nextStep(void);                    //take next step in simulation
void gridColor(int i, int j);           //change grid color at each node
void myDisplay(void);                   //GLUT display function
void myReshape(int w, int h);           //GLUT reshape function
void myKeyboard(unsigned char key, int xx, int yy); //GLUT keyboard function
void mySpecialKey(int key, int xx, int yy);         //GLUT special keys function
void exitProgram(void);                 //free memory and exit program
void clRhs(void);                       //call kernel for calculating RHS
int setupCL(const char *filename);      //setup opencl framework
char* loadProgramSource(const char *filename); //load kernel source code from file

int main (int argc, char *argv[])
{
    printf("**************************** HYBRID 2D ****************************\n");
    readAndCheckInput();
    initialConditions();
    printOutput();
    printf("Running...");
    fflush(stdout);
    
    //visualization setup for GLUT
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB|GLUT_DEPTH);
    glutInitWindowSize(winx, winy);
    glutCreateWindow("HYBRID2D");
    glutDisplayFunc(myDisplay);
    glutIdleFunc(myDisplay);
    glutReshapeFunc(myReshape);
    glutKeyboardFunc(myKeyboard);
    glutSpecialFunc(mySpecialKey);
    glMatrixMode(GL_MODELVIEW);
    glClearColor(0.0f,0.0f,0.0f,1.0f);
    glEnable(GL_DEPTH_TEST);
    glutMainLoop();  //run simulation
    
    return 0;
}

char* loadProgramSource(const char *filename)
{
    struct stat statbuf;
    FILE *fh;
    char *source;
    
    fh = fopen(filename, "r");
    if (fh == 0)
        return 0;
    
    stat(filename, &statbuf);
    source = (char *) malloc(statbuf.st_size + 1);
    fread(source, statbuf.st_size, 1, fh);
    source[statbuf.st_size] = '\0';
    fclose(fh);
    return source;
}

int setupCL(const char *filename)
{
    //connect to a compute device
    err = clGetDeviceIDs(NULL,CL_DEVICE_TYPE_CPU, 1, &devices, NULL);
    //err = clGetDeviceIDs(NULL,CL_DEVICE_TYPE_GPU, 1, &devices, NULL);
    
    //get info about device
    size_t returned_size = 0;
    cl_char vendor_name[1024] = {0};
    cl_char device_name[1024] = {0};
    err = clGetDeviceInfo(devices, CL_DEVICE_VENDOR, sizeof(vendor_name), vendor_name, &returned_size);
    err|= clGetDeviceInfo(devices, CL_DEVICE_NAME, sizeof(device_name), device_name, &returned_size);
    printf("Connecting to %s %s...\n", vendor_name, device_name);
    
    //read the program
    printf("Loading program '%s'\n\n", filename);
    char *program_source = loadProgramSource(filename);
    
    //create the context and command queue
    context = clCreateContext(0, 1, &devices, NULL, NULL, &err);
    cmd_queue = clCreateCommandQueue(context, devices, 0, NULL);
    
    //create program from .cl file
    program = clCreateProgramWithSource(context,1, (const char**)&program_source, NULL, &err);
    
    //build the kernel program
    err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    printf("Error: %d\n", err);
    char build[2048];
    clGetProgramBuildInfo(program, devices, CL_PROGRAM_BUILD_LOG, 2048, build, NULL);
    printf("Build Log:\n%s\n",build);
    
    //create the kernel
    kernel = clCreateKernel(program, "rhs", &err);
    
    //get size of buffers
    size_t w_n_buffer_size = sizeof(float)*KSIZE*X*Y;
    size_t K_buffer_size = sizeof(float)*KSIZE*X*Y;
    size_t p_buffer_size = sizeof(float)*X*Y;
    size_t H0_buffer_size = sizeof(float)*X*Y;
    
    //create buffers
    w_n_mem = clCreateBuffer(context, CL_MEM_READ_ONLY,  w_n_buffer_size, NULL, NULL);
    p_mem   = clCreateBuffer(context, CL_MEM_READ_ONLY,  p_buffer_size,   NULL, NULL);
    H0_mem  = clCreateBuffer(context, CL_MEM_READ_ONLY,  H0_buffer_size,  NULL, NULL);
    K_mem   = clCreateBuffer(context, CL_MEM_READ_WRITE, K_buffer_size,   NULL, NULL);
    
    //set kernel arguments
    err  = clSetKernelArg(kernel,  0, sizeof(cl_mem), &w_n_mem);
    err |= clSetKernelArg(kernel,  1, sizeof(cl_mem), &K_mem);
    err |= clSetKernelArg(kernel,  2, sizeof(cl_mem), &p_mem);
    err |= clSetKernelArg(kernel,  3, sizeof(cl_mem), &H0_mem);
    err |= clSetKernelArg(kernel,  4, sizeof(int),    &Y);
    err |= clSetKernelArg(kernel,  5, sizeof(int),    &X);
    err |= clSetKernelArg(kernel,  6, sizeof(float),  &dy);
    err |= clSetKernelArg(kernel,  7, sizeof(float),  &dx);
    err |= clSetKernelArg(kernel,  8, sizeof(float),  &alpha);
    
    return CL_SUCCESS;
}

void clRhs(void)
{
    //get size of buffers
    size_t w_n_buffer_size = sizeof(float)*KSIZE*X*Y;
    size_t K_buffer_size = sizeof(float)*KSIZE*X*Y;
    size_t p_buffer_size = sizeof(float)*X*Y;
    size_t H0_buffer_size = sizeof(float)*X*Y;
    
    //set work-item dimensions
    size_t global_work_size[2];
    global_work_size[0] = X; // Number of computation to perform.
    global_work_size[1] = Y; // Number of computation to perform.
    
    //queue memory to be written to the device
    err = clEnqueueWriteBuffer(cmd_queue, w_n_mem, CL_TRUE, 0, w_n_buffer_size, (void*)w_n, 0, NULL, NULL);
    err = clEnqueueWriteBuffer(cmd_queue, p_mem,   CL_TRUE, 0, p_buffer_size,   (void*)p,   0, NULL, NULL);
    err = clEnqueueWriteBuffer(cmd_queue, H0_mem,  CL_TRUE, 0, H0_buffer_size,  (void*)H0,  0, NULL, NULL);
    clFinish(cmd_queue);
    
    //run calculation
    err = clEnqueueNDRangeKernel(cmd_queue, kernel, 2, NULL, global_work_size, NULL, 0, NULL, NULL);
    clFinish(cmd_queue);
    
    //get results from device
    err = clEnqueueReadBuffer(cmd_queue, K_mem, CL_TRUE, 0, K_buffer_size, K, 0, NULL, NULL);
    clFinish(cmd_queue);
}

void exitProgram(void)
{
    printf("\rRunning: Done.                               \n");
    freeMemory();
    exit(0);
}

void gridColor(int i, int j)
{
    float rr, gg, bb;
    
    rr = (float)(fabs(b[(j*X)+i])/bmax);
    gg = 1.0f;
    bb = (float)(b[(j*X)+i]/bmax);
    
    glColor3f(rr, gg, bb);
}

void myDisplay(void)
{
    if (!paused) {
        if (step < time_steps) {
            step++;
            nextStep();
        } else {
            exitProgram();
        }
    }
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    gluLookAt((X-1.0)/2.0, -Y/8.0, (X+Y)/2.6, (X-1.0)/grid_look_x, Y/grid_look_y, 0.0, 0.0, 1.0, 0.0);
    bmax = 0;
    for (int j = 1; j < Y-1; j++) {
        for (int i = 1; i < X-1; i++) {
            b[(j*X)+i] = grid_amp*(p[(j*X)+i]-p0);
            if (b[(j*X)+i] > bmax)
                bmax = b[(j*X)+i];    //set global maximum for grid coloring
        }
    }
    //draw grid
    glPushMatrix();
    glTranslatef((float)X/2, (float)Y/2,0);
    glRotatef((float)grid_rot_x, 1.0f, 0.0f, 0.0f);
    glRotatef((float)grid_rot_z, 0.0f, 0.0f, 1.0f);
    glTranslatef(-(float)X/2, -(float)Y/2, 0);
    for (int j = 0; j < Y-1; j++) {
        for (int i = 0; i < X-1; i++) {
            glBegin(GL_LINES);
            gridColor(i, j);
            glVertex3f((float)i,   (float)j, (float)b[(j*X)+i]);
            gridColor(i+1, j);
            glVertex3f((float)i+1, (float)j, (float)b[(j*X)+i+1]);
            glEnd();
            
            glBegin(GL_LINES);
            gridColor(i, j);
            glVertex3f((float)i, (float)j,   (float)b[(j*X)+i]);
            gridColor(i, j+1);
            glVertex3f((float)i, (float)j+1, (float)b[((j+1)*X)+i]);
            glEnd();
        }
        glBegin(GL_LINES);
        gridColor(X-1, j);
        glVertex3f((float)X-1, (float)j,   (float)b[(j*X)+X-1]);
        gridColor(X-1, j+1);
        glVertex3f((float)X-1, (float)j+1, (float)b[((j+1)*X)+X-1]);
        glEnd();
    }
    for (int i = 0; i < X-1; i++) {
        glBegin(GL_LINES);
        gridColor(i, Y-1);
        glVertex3f((float)i,   (float)Y-1, (float)b[((Y-1)*X)+i]);
        gridColor(i+1, Y-1);
        glVertex3f((float)i+1, (float)Y-1, (float)b[((Y-1)*X)+i+1]);
        glEnd();
    }
    glPopMatrix();
    glutSwapBuffers();
}

void myReshape(int w, int h)
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(grid_fov, 1.0*w/h, 1, Y*2);
    glViewport(0, 0, w, h);
    glMatrixMode(GL_MODELVIEW);
}

void myKeyboard(unsigned char key, int xx, int yy)
{
    switch (key) {
        case 13 :
            if (paused == 1) {paused = 0;}
            else if (paused == 0) {paused = 1;}
            break;
        case 'q' :
            if (pabsorb == 1) {pabsorb = 0;}
            else if (pabsorb == 0) {pabsorb = 1;}
            break;
        case ']' :
            grid_amp *= 2;
            break;
        case '[' :
            grid_amp /= 2;
            break;
        case '.' :
            grid_fov -= 1;
            myReshape(winx, winy);
            break;
        case ',' :
            grid_fov += 1;
            myReshape(winx, winy);
            break;
        case 'w' :
            grid_look_y -= 0.05;
            break;
        case 's' :
            grid_look_y += 0.05;
            break;
        case 'a' :
            grid_look_x += 0.05;
            break;
        case 'd' :
            grid_look_x -= 0.05;
            break;
    }
}

void mySpecialKey(int key, int xx, int yy)
{
    switch (key) {
        case GLUT_KEY_UP :
            grid_rot_x -= 2;
            break;
        case GLUT_KEY_DOWN :
            grid_rot_x += 2;
            break;
        case GLUT_KEY_LEFT :
            grid_rot_z += 2;
            break;
        case GLUT_KEY_RIGHT :
            grid_rot_z -= 2;
            break;
    }
}

void initialConditions(void)
{
    lambda = c/frequency_max;    //minimum wavelength to account for
    dx = lambda/ppw;             //delta x
    dy = lambda/ppw;             //delta y
    dt = cfl*(dx/c);             //delta t
    x_length = dx*X;             //length of x domain in meters
    y_length = dy*Y;             //length of y domain in meters
    
    b = (float *)malloc(Y*X*sizeof(float)); if (b == NULL) mfail("b");
    p = (float*)malloc(X*Y*sizeof(float)); if (p == NULL) mfail("p");
    T = (float*)malloc(X*Y*sizeof(float)); if (T == NULL) mfail("T");
    H0 = (float*)malloc(X*Y*sizeof(float)); if (H0 == NULL) mfail("H0");
    w = (float*)malloc(KSIZE*X*Y*sizeof(float)); if (w == NULL) mfail("w");
    w_n = (float*)malloc(KSIZE*X*Y*sizeof(float)); if (w_n == NULL) mfail("w_n");
    K = (float*)malloc(KSIZE*X*Y*sizeof(float)); if (K == NULL) mfail("K");
    
    //initialize grids to ambient conditions
    #pragma omp parallel for
    for (int j = 0; j < Y; j++) {
        for (int i = 0; i < X; i++) {
            b[(j*X)+i] = 0;
            p[(j*X)+i] = p0;
            T[(j*X)+i] = T0-(0.01*j*dy);
            H0[(j*X)+i] = 0;
        }
    }
    
    //set all vectors for equations to 0; if using calloc this can go away
    #pragma omp parallel for
    for (int j = 0; j < Y; j++) {
        for (int k = 0; k < KSIZE; k++) {
            for (int i = 0; i < X; i++) {
                w[(k*Y*X)+(j*X)+i] = w_n[(k*Y*X)+(j*X)+i] = K[(k*Y*X)+(j*X)+i] = 0;
            }
        }
    }
    
    //calculate relaxation frequencies and times
    const float zeta = 10.79586*(1-(273.16/T0))-5.02808*log10(T0/273.16)
        +(1.50474e-4)*(1-pow(10,(-8.29692*((T0/273.16)-1))))
        -(4.2873e-4)*(1-pow(10,(-4.76955*((273.16/T0)-1))))-2.2195983;
    const float p_vp = p_ref*pow(10,zeta);
    const float h = (pow(10.0,-2)*RH*p_vp*T0)/p0;
    const float eta = 4.17*(pow(T_ref/T0,1.0/3.0)-1);
    const float f_n = (p0/p_ref)*(24+(4.04e6)*h*((0.02+100*h)/(0.391+100*h)));
    const float f_o = (p0/p_ref)*(pow((T_ref/T0),0.5)*(9+2.8e4*h*exp(-eta)));
    tau_n = 1/(2*PI*f_n);
    tau_o = 1/(2*PI*f_o);
    
    //initialize w to ambient values
    for (int j = 0; j < Y; j++) {
        for (int i = 0; i < X; i++) {
            w[(0*Y*X)+(j*X)+i] = rho0;
            w[(1*Y*X)+(j*X)+i] = rho0*0;
            w[(2*Y*X)+(j*X)+i] = rho0*0;
            w[(3*Y*X)+(j*X)+i] = rho0*s0;
            //w[(4*Y*X)+(j*X)+i] = rho0*T0;
            //w[(5*Y*X)+(j*X)+i] = rho0*T0;
        }
    }
    
    //set w_n=w
    #pragma omp parallel for
    for (int j = 0; j < Y; j++) {
        for (int k = 0; k < KSIZE; k++) {
            for (int i = 0; i < X; i++) {
                w_n[(k*Y*X)+(j*X)+i] = w[(k*Y*X)+(j*X)+i];
            }
        }
    }
    
    updateVariables();
    
    setupCL("kernel.cl");
}

void nextStep(void)
{
    //if source is a pulse, add it only at step 1, if driven, always add it in
    if (source_mode == 0) {
        if (step == 1) {
            addSource();
        }
    } else if (source_mode == 1) {
        addSource();
    }
    
    //perform runge kutta
    rungeKutta();
    
    //if source is a pulse, reset H density term to 0
    if (source_mode == 0) {
        for (int j = 0; j < Y; j++) {
            for (int i = 0; i < X; i++) {
                H0[(j*X)+i] = 0;
            }
        }
    }
}

void updateVariables(void)
{
    //update p and T
    #pragma omp parallel for
    for (int j = 0; j < Y; j++) {
        for (int i = 0; i < X; i++) {
            p[(j*X)+i] = c*c*((w_n[(0*Y*X)+(j*X)+i]-rho0)+((gamma1-1)/(2*rho0))*SQUARE(w_n[(0*Y*X)+(j*X)+i]-rho0)
                +(w_n[(0*Y*X)+(j*X)+i]*(T[(j*X)+i]/T0)/c_p)*((w_n[(3*Y*X)+(j*X)+i]/w_n[(0*Y*X)+(j*X)+i])-s0))+p0;
            T[(j*X)+i] = (T[(j*X)+i]/c_p)*((w_n[(3*Y*X)+(j*X)+i]/w_n[(0*Y*X)+(j*X)+i])-s0)+((T[(j*X)+i]/T0)/(w_n[(0*Y*X)+(j*X)+i]*c_p))
                *(p[(j*X)+i]-p0)+T0;
        }
    }
}

void rungeKutta(void)
{
    //stage 1
    clRhs();
    #pragma omp parallel for
    for (int j = 0; j < Y; j++) {
        for (int k = 0; k < KSIZE; k++) {
            for (int i = 0; i < X; i++) {
                w_n[(k*Y*X)+(j*X)+i] = w[(k*Y*X)+(j*X)+i]+dt*K[(k*Y*X)+(j*X)+i];
            }
        }
    }
    updateVariables();
    
    //stage 2
    clRhs();
    #pragma omp parallel for
    for (int j = 0; j < Y; j++) {
        for (int k = 0; k < KSIZE; k++) {
            for (int i = 0; i < X; i++) {
                w_n[(k*Y*X)+(j*X)+i] = 0.5*w[(k*Y*X)+(j*X)+i]+0.5*w_n[(k*Y*X)+(j*X)+i]+0.5*dt*K[(k*Y*X)+(j*X)+i];
                w[(k*Y*X)+(j*X)+i] = w_n[(k*Y*X)+(j*X)+i];
            }
        }
    }
    updateVariables();
    
    //apply boundary conditions
    bc();
}

void addSource(void)
{
    if (source_type == 1) {
        if (source_mode == 1) {
            for (int j = 0; j < Y; j++) {
                for (int i = 0; i < X; i++) {
                    //gaussian distributed smooth point source
                    H0[(j*X)+i] = source_amp*exp(-log(2)/SQUARE(alph*dx)
                        *(SQUARE(i+1-src_x_pos)+SQUARE(j+1-src_y_pos)))
                        *sin(omega*step*dt);
                }
            }
        } else if (source_mode == 0) {
            for (int j = 0; j < Y; j++) {
                for (int i = 0; i < X; i++) {
                    //gaussian distributed smooth point source
                    H0[(j*X)+i] = source_amp*exp(-log(2)/SQUARE(alph*dx)
                        *(SQUARE(i+1-src_x_pos)+SQUARE(j+1-src_y_pos)));
                }
            }
        }
    } else if (source_type == 2) {
        if (source_mode == 1) {
            for (int j = 0; j < Y; j++) {
                for (int i = 0; i < X; i++) {
                    //sine wave smooth source
                    H0[(j*X)+i] = source_amp*exp(-alph*dx
                        *SQUARE((i+1-src_x_pos)-sine_src_amp*sin(sine_src_frq*(j+1-src_y_pos))))
                        *sin(omega*step*dt);
                }
            }
        } else if (source_mode == 0) {
            for (int j = 0; j < Y; j++) {
                for (int i = 0; i < X; i++) {
                    //sine wave smooth source
                    H0[(j*X)+i] = source_amp*exp(-alph*dx
                        *SQUARE((i+1-src_x_pos)-sine_src_amp*sin(sine_src_frq*(j+1-src_y_pos))));
                }
            }
        }
    }
}

void bc(void)
{
    const int numpoints = 10;        //number of points from boundary to start absorbing from
    const float x_p = 0;             //numpoints*dx;
    const float alx = 10*ppw*dx;     //strength of absorption
    const float y_p = 0;             //numpoints*dy;
    const float aly = 10*ppw*dy;     //strength of absorption
    
    if (bc_left==1) {
        //set column boundaries equal to next column
        for (int j = 0; j < Y; j++) {
            for (int k = 0; k < KSIZE; k++) {
                w[(k*Y*X)+(j*X)+0] = w[(k*Y*X)+(j*X)+1]=w[(k*Y*X)+(j*X)+2]=w[(k*Y*X)+(j*X)+3];
            }
        }
        
        //reflect off left side
        for (int j=0; j<Y; j++) {
            w[(1*Y*X)+(j*X)+2]=-w[(1*Y*X)+(j*X)+2];
        }
    } else if (bc_left == 2) { //absorb
        for (int j = 0; j < Y; j++) {
            for (int i = 0; i < numpoints; i++) {
                float tmp = (-exp(-log(2.0)/SQUARE(alx)*SQUARE(i-x_p))+1);
                w[(0*Y*X)+(j*X)+i] = (w[(0*Y*X)+(j*X)+i]-rho0)*tmp+rho0;
                w[(1*Y*X)+(j*X)+i] = w[(1*Y*X)+(j*X)+i]*tmp;
                w[(2*Y*X)+(j*X)+i] = w[(2*Y*X)+(j*X)+i]*tmp;
                w[(3*Y*X)+(j*X)+i] = w[(3*Y*X)+(j*X)+i]*tmp;
                //w[(4*Y*X)+(j*X)+i] = (w[(4*Y*X)+(j*X)+i]-T0*rho0)*tmp+T0*rho0;
                //w[(5*Y*X)+(j*X)+i] = (w[(5*Y*X)+(j*X)+i]-T0*rho0)*tmp+T0*rho0;
                p[(j*X)+i] = (p[(j*X)+i]-p0)*tmp+p0;
                T[(j*X)+i] = (T[(j*X)+i]-T0)*tmp+T0;
            }
        }
    } else if (bc_left == 3) {
        for (int j = 0; j < Y; j++) {
            for (int k = 0; k < KSIZE; k++) {
                for (int i = 0; i < 10; i++) {
                    w[(k*Y*X)+(j*X)+X-10+i] = w[(k*Y*X)+(j*X)+10+i];
                    w[(k*Y*X)+(j*X)+i] = w[(k*Y*X)+(j*X)+X-20+i];
                }
            }
        }
        for (int j = 0; j < Y; j++) {
            for (int i = 0; i < 10; i++) {
                p[(j*X)+X-10+i] = p[(j*X)+10+i];
                p[(j*X)+i] = p[(j*X)+X-20+i];
                T[(j*X)+X-10+i] = T[(j*X)+10+i];
                T[(j*X)+i] = T[(j*X)+X-20+i];
            }
        }
    }
    
    if (bc_right == 1) {
        //set column boundaries equal to next column
        for (int j = 0; j < Y; j++) {
            for (int k = 0; k < KSIZE; k++) {
                w[(k*Y*X)+(j*X)+X-1] = w[(k*Y*X)+(j*X)+X-2]=w[(k*Y*X)+(j*X)+X-3]=w[(k*Y*X)+(j*X)+X-4];
            }
        }
        
        //reflect off right side
        for (int j = 0; j < Y; j++) {
            w[(1*Y*X)+(j*X)+X-3] = -w[(1*Y*X)+(j*X)+X-3];
        }
    } else if (bc_right == 2) {
        for (int j = 0; j < Y; j++) {
            for (int i = X-numpoints; i < X; i++) {
                float tmp = (-exp(-log(2.0)/SQUARE(alx)*SQUARE(X-i-1-x_p))+1);
                w[(0*Y*X)+(j*X)+i] = (w[(0*Y*X)+(j*X)+i]-rho0)*tmp+rho0;
                w[(1*Y*X)+(j*X)+i] = w[(1*Y*X)+(j*X)+i]*tmp;
                w[(2*Y*X)+(j*X)+i] = w[(2*Y*X)+(j*X)+i]*tmp;
                w[(3*Y*X)+(j*X)+i] = w[(3*Y*X)+(j*X)+i]*tmp;
                //w[(4*Y*X)+(j*X)+i] = (w[(4*Y*X)+(j*X)+i]-T0*rho0)*tmp+T0*rho0;
                //w[(5*Y*X)+(j*X)+i] = (w[(5*Y*X)+(j*X)+i]-T0*rho0)*tmp+T0*rho0;
                p[(j*X)+i] = (p[(j*X)+i]-p0)*tmp+p0;
                T[(j*X)+i] = (T[(j*X)+i]-T0)*tmp+T0;
            }
        }
    } else if (bc_right == 3) {
        for (int j = 0; j < Y; j++) {
            for (int k = 0; k < KSIZE; k++) {
                for (int i = 0; i < 10; i++) {
                    w[(k*Y*X)+(j*X)+i] = w[(k*Y*X)+(j*X)+X-20+i];
                    w[(k*Y*X)+(j*X)+X-10+i] = w[(k*Y*X)+(j*X)+10+i];
                }
            }
        }
        for (int j = 0; j < Y; j++) {
            for (int i = 0; i < 10; i++) {
                p[(j*X)+i] = p[(j*X)+X-20+i];
                p[(j*X)+X-10+i] = p[(j*X)+10+i];
                T[(j*X)+i] = T[(j*X)+X-20+i];
                T[(j*X)+X-10+i] = T[(j*X)+10+i];
            }
        }
        if (pabsorb) {
            for (int j = 0; j < Y; j++) {
                for (int i = 0; i < numpoints; i++) {
                    float tmp = (-exp(-log(2.0)/SQUARE(alx)*SQUARE(i-x_p))+1);
                    w[(0*Y*X)+(j*X)+i+30] = (w[(0*Y*X)+(j*X)+i+30]-rho0)*tmp+rho0;
                    w[(1*Y*X)+(j*X)+i+30] = w[(1*Y*X)+(j*X)+i+30]*tmp;
                    w[(2*Y*X)+(j*X)+i+30] = w[(2*Y*X)+(j*X)+i+30]*tmp;
                    w[(3*Y*X)+(j*X)+i+30] = w[(3*Y*X)+(j*X)+i+30]*tmp;
                    //w[(4*Y*X)+(j*X)+i+30] = (w[(4*Y*X)+(j*X)+i+30]-T0*rho0)*tmp+T0*rho0;
                    //w[(5*Y*X)+(j*X)+i+30] = (w[(5*Y*X)+(j*X)+i+30]-T0*rho0)*tmp+T0*rho0;
                    p[(j*X)+i+30] = (p[(j*X)+i+30]-p0)*tmp+p0;
                    T[(j*X)+i+30] = (T[(j*X)+i+30]-T0)*tmp+T0;
                }
            }
        }
    }
    
    if (bc_bottom == 1) {
        //set row boundaries equal to next row
        for (int i = 0; i < X; i++) {
            for (int k = 0; k < KSIZE; k++) {
                w[(k*Y*X)+(0*X)+i] = w[(k*Y*X)+(1*X)+i]=w[(k*Y*X)+(2*X)+i]=w[(k*Y*X)+(3*X)+i];
            }
        }
        
        //reflect off bottom side
        for (int i = 0; i < X; i++) {
            w[(2*Y*X)+(2*X)+i] = -w[(2*Y*X)+(2*X)+i];
        }
    } else if (bc_bottom == 2) {
        for (int j = 0; j < numpoints; j++) {
            for (int i = 0; i < X; i++) {
                float tmp = (-exp(-log(2.0)/SQUARE(aly)*SQUARE(j-y_p))+1);
                w[(0*Y*X)+(j*X)+i] = (w[(0*Y*X)+(j*X)+i]-rho0)*tmp+rho0;
                w[(1*Y*X)+(j*X)+i] = w[(1*Y*X)+(j*X)+i]*tmp;
                w[(2*Y*X)+(j*X)+i] = w[(2*Y*X)+(j*X)+i]*tmp;
                w[(3*Y*X)+(j*X)+i] = w[(3*Y*X)+(j*X)+i]*tmp;
                //w[(4*Y*X)+(j*X)+i] = (w[(4*Y*X)+(j*X)+i]-T0*rho0)*tmp+T0*rho0;
                //w[(5*Y*X)+(j*X)+i] = (w[(5*Y*X)+(j*X)+i]-T0*rho0)*tmp+T0*rho0;
                p[(j*X)+i] = (p[(j*X)+i]-p0)*tmp+p0;
                T[(j*X)+i] = (T[(j*X)+i]-T0)*tmp+T0;
            }
        }
    } else if (bc_bottom == 3) {
        for (int k = 0; k < KSIZE; k++) {
            for (int j = 0; j < 10; j++) {
                for (int i = 0; i < X; i++) {
                    w[(k*Y*X)+((Y-10+j)*X)+i] = w[(k*Y*X)+((10+j)*X)+i];
                    w[(k*Y*X)+(j*X)+i] = w[(k*Y*X)+((Y-20+j)*X)+i];
                }
            }
        }
        for (int j = 0; j < 10; j++) {
            for (int i = 0; i < X; i++) {
                p[((Y-10+j)*X)+i] = p[((10+j)*X)+i];
                p[(j*X)+i] = p[((Y-20+j)*X)+i];
                T[((Y-10+j)*X)+i] = T[((10+j)*X)+i];
                T[(j*X)+i] = T[((Y-20+j)*X)+i];
            }
        }
    }
    
    if (bc_top == 1) {
        //set row boundaries equal to next row
        for (int k = 0; k < KSIZE; k++) {
            for (int i = 0; i < X; i++) {
                w[(k*Y*X)+((Y-1)*X)+i] = w[(k*Y*X)+((Y-2)*X)+i]=w[(k*Y*X)+((Y-3)*X)+i]=w[(k*Y*X)+((Y-4)*X)+i];
            }
        }
        
        //reflect off top side
        for (int i = 0; i < X; i++) {
            w[(2*Y*X)+((Y-3)*X)+i] = -w[(2*Y*X)+((Y-3)*X)+i];
        }
    } else if (bc_top == 2) {
        for (int j = Y-numpoints; j < Y; j++) {
            for (int i = 0; i < X; i++) {
                float tmp = (-exp(-log(2.0)/SQUARE(aly)*SQUARE(Y-j-1-y_p))+1);
                w[(0*Y*X)+(j*X)+i] = (w[(0*Y*X)+(j*X)+i]-rho0)*tmp+rho0;
                w[(1*Y*X)+(j*X)+i] = w[(1*Y*X)+(j*X)+i]*tmp;
                w[(2*Y*X)+(j*X)+i] = w[(2*Y*X)+(j*X)+i]*tmp;
                w[(3*Y*X)+(j*X)+i] = w[(3*Y*X)+(j*X)+i]*tmp;
                //w[(4*Y*X)+(j*X)+i] = (w[(4*Y*X)+(j*X)+i]-T0*rho0)*tmp+T0*rho0;
                //w[(5*Y*X)+(j*X)+i] = (w[(5*Y*X)+(j*X)+i]-T0*rho0)*tmp+T0*rho0;
                p[(j*X)+i] = (p[(j*X)+i]-p0)*tmp+p0;
                T[(j*X)+i] = (T[(j*X)+i]-T0)*tmp+T0;
            }
        }
    } else if (bc_top == 3) {
        for (int k = 0; k < KSIZE; k++) {
            for (int j = 0; j < 10; j++) {
                for (int i = 0; i < X; i++) {
                    w[(k*Y*X)+(j*X)+i] = w[(k*Y*X)+((Y-20+j)*X)+i];
                    w[(k*Y*X)+((Y-10+j)*X)+i] = w[(k*Y*X)+((10+j)*X)+i];
                }
            }
        }
        for (int j = 0; j < 10; j++) {
            for (int i = 0; i < X; i++) {
                p[(j*X)+i] = p[((Y-20+j)*X)+i];
                p[((Y-10+j)*X)+i] = p[((10+j)*X)+i];
                T[(j*X)+i] = T[((Y-20+j)*X)+i];
                T[((Y-10+j)*X)+i] = T[((10+j)*X)+i];
            }
        }
    }
}

void freeMemory(void)
{
    free(p);
    free(T);
    free(H0);
    free(b);
    free(w);
    free(w_n);
    free(K);
    
    // release kernel, program, and memory objects
    clReleaseKernel(kernel);
    clReleaseProgram(program);
    clReleaseCommandQueue(cmd_queue);
    clReleaseContext(context);
    clReleaseMemObject(w_n_mem);
    clReleaseMemObject(K_mem);
    clReleaseMemObject(p_mem);
    clReleaseMemObject(H0_mem);
}

void fail(char *message)
{
    printf("%s. Goodbye.", message);
    exit(1);
}

void mfail(char *message)
{
    printf("Error allocating memory for %s. Goodbye.\n", message);
    exit(1);
}

void printOutput(void)
{
    printf("Maximum Frequency: %d Hz, Minimum Wavelength: %g m\n", frequency_max, lambda);
    printf("Grid Cells X: %d, Domain Length X: %g m, dx: %g m\n", X, x_length, dx);
    printf("Grid Cells Y: %d, Domain Length Y: %g m, dy: %g m\n", Y, y_length, dy);
    printf("Time Steps: %d, Total Time: %g s, dt: %g s\n", time_steps, time_steps*dt,dt);
    printf("CFL: %g, PPW: %d, Splitting alpha: %g, gthres: %g\n", cfl, ppw, alpha, gthres);
    printf("Source alpha: %g, Source omega: %g, Source Pos: (%d,%d)\n", alph, omega, src_x_pos, src_y_pos);
    printf("Source amplitude: %g, Source mode: %d\n", source_amp, source_mode);
}

void readAndCheckInput(void)
{
    if ((winx < 10) || (winy < 10))
        fail("That's a pretty small window to create.");
    if ((winx > 3000) || (winy > 3000))
        fail("That's a pretty large window to create.");
    if ((cfl <= 0) || (cfl >= 0.5))
        fail("CFL is out of bounds.");
    if (X < 21)
        fail("X is too small or negative.");
    if (Y < 21)
        fail("Y is too small or negative.");
    if (time_steps < 1)
        fail("Number of time steps is too small or negative.");
    if ((bc_left < 0) || (bc_left > 3))
        fail("BC left value is not 0,1,2 or 3.");
    if ((bc_right < 0) || (bc_right > 3))
        fail("BC right value is not 0,1,2 or 3.");
    if ((bc_bottom < 0) || (bc_bottom > 3))
        fail("BC bottom value is not 0,1,2 or 3.");
    if ((bc_top < 0) || (bc_top > 3))
        fail("BC top value is not 0,1,2 or 3.");
    if (ppw < 10)
        fail("PPW value is too small.");
    if (frequency_max < 0)
        fail("Maximum frequency is negative.");
    if ((src_x_pos < 0) || (src_x_pos > X))
        fail("Source x position is outside the domain.");
    if ((src_y_pos < 0) || (src_y_pos > Y))
        fail("Source y position is outside the domain.");
    if ((source_mode < 0) || (source_mode > 1))
        fail("Source mode is neither 0 nor 1.");
    if (alph < 1)
        fail("Source width is 0 or negative.");
    if (alpha < 300)
        fail("Flux splitting alpha is too small.");
    if ((gthres < 0) || (gthres > 1))
        fail("Threshold should be between 0 and 1.");
    if ((source_type != 1) && (source_type != 2) && (source_type != 3))
        fail("Source type should be 1, 2 or 3.");
    if ((X > 1000) || (Y > 1000))
        fail("Viewing too large of a domain.");
    
    printf("Successfully checked input parameters.\n"); fflush(stdout);
}
