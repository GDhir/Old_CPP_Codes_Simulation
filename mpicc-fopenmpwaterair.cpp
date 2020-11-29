#include<iostream>
#include<cmath>
#include<fstream>
#include<sstream>
#include<string>
# define ddd 10
using namespace std;
// This structure represents the x and y components of each vector
struct vectorc
{
    double x;
    double y;
};

/* The following structure consists of properties defined within each cell for each phase. In a 2D rectangular grid, we have four
surfaces bounding the cell. Hence the quantities defined at interface or boundary will be represented as an array of four
values inside the cell.*/
struct cell
{
    vectorc position[4]; // Coordinates of four vertices.
    vectorc center;        // Coordinates of the center.
    vectorc area[4];       // Vector representing the four surface areas.
    double volumet[4];     // We divide the cell in four areas, with each one consisting of a triangle bounded by two successive vertices and the center of cell.
    vectorc nn[4];         // Unit vector Normal to each area
    vectorc pp[4];         // Unit Vector tangent to the four sides of cell.
    vectorc side[4];       // Vector representing four sides of the cell
    double volume;         // Volume of the cell
    double r[4];
    double Cp, pinf, gamma;  // Constants used in equation of state
    double alpha;            // Volume fraction of each phase
    double alphal[4];         //  Reconstructed Volume fraction on the left
    double alphar[4];         // Reconstructed Volume fraction on right
    double p;                 // Pressure
    double pl[4];             // Reconstructed pressure on the left
    double pr[4];             //  Reconstructed pressure on the right
    double pint;              //  Pressure at the interface between two phases inside the cell
    double E;                  // Total Energy per unit mass
    double e;                  // Internal Energy per unit mass
    double el[4];              //  Reconstructed internal energy on the left
    double er[4];              //  Reconstructed internal energy on the right
    double T;                   // Temperature
    double rho;                 // density
    double rhol[4];             //  Reconstructed density on the left
    double rhor[4];             //  Reconstructed density on the right
    vectorc vcart;              //  Velocity in each cell in cartesian coordinates
    vectorc vlocal[4];          // Velocity in each cell in local coordinates where local coordinates refer to coordinate system of each boundary area.
    vectorc vllocal[4];         //  Reconstructed velocity in local coordinates on the left
    vectorc vrlocal[4];         //  Reconstructed velocity in local coordinates on the right
    vectorc vavg;               //  Average velocity of both phases in each cell.
    double Tavg;                // Average Temperature of both phases in each cell.
    double rhoavg;              // Average density of both phases
    double a;                   // Speed of sound
    double al[4];               //  Reconstructed speed of sound on the left
    double ar[4];               //  Reconstructed speed of sound on the right
    double ic_ak[4];            //  Speed of sound for each phase on the boundary surfaces of each cell
    double ic_a[4];
    double H;                    //  Enthalpy
    double Hl[4],Hr[4],Tl[4],Tr[4];    // Reconstructed values on the left and right
    double deltaeff[4];                // Represents phase transition of each phase across each bounding surface of the cell
    double U[4];                       // Conservative variables [alpha*rho,alpha*rho*u,alpha*rho*v,alpha*rho*E]
    double W[4];
    vectorc ifc[4];                    //Coordinates of Center of each bounding area
    vectorc delrho;                    //  Density gradient in each cell
    double widx;
    double widy;
};

//Following functions represent vector operations
double crossprod(vectorc a,vectorc b)    // Returns Crossproduct of two vectors
{
    double c;
    c=a.x*b.y-a.y*b.x;
    return c;
}
vectorc subtract(vectorc a,vectorc b)   //  Subtracts two vectors
{
    vectorc c;
    c.x=a.x-b.x;
    c.y=a.y-b.y;
    return c;
}
vectorc add(vectorc a,vectorc b)          //  Adds two vectors
{
    vectorc c;
    c.x=a.x+b.x;
    c.y=a.y+b.y;
    return c;
}
double mod1(vectorc a)       //  Modulus of each vector
{
    double c;
    c=sqrt(pow(a.x,2)+pow(a.y,2));
    return c;
}
double dot(vectorc a,vectorc b)     //  dot product of two vectors
{
    double c;
    c=a.x*b.x+a.y*b.y;
    return c;
}

// Following function computes geometric variables associated with each cell
void geometry(cell***b,int nx,int ny)
{
    vectorc d;
    int t;
    for(int l=0;l<2;l++)
    {
        for(int i=0;i<(nx-1);i++)
        {
            for(int j=0;j<(ny-1);j++)
            {
                b[i][j][l].center.x=0;
                b[i][j][l].center.y=0;
                for(int h=0;h<4;h++)
                {
                    b[i][j][l].center=add(b[i][j][l].position[h],b[i][j][l].center);
                }
                b[i][j][l].center.x=b[i][j][l].center.x/4;
                b[i][j][l].center.y=b[i][j][l].center.y/4;
                b[i][j][l].volume=0;
                // Evaluating volume of each cell
                for(int h=0;h<4;h++)
                {
                    if(h==3)
                        t=0;
                    else
                        t=h+1;
                    b[i][j][l].volumet[h]=0.5*crossprod(subtract(b[i][j][l].position[h],b[i][j][l].center),subtract(b[i][j][l].position[t],b[i][j][l].center));
                    b[i][j][l].volume=b[i][j][l].volume+b[i][j][l].volumet[h];
                    //cout<<area<<"\n";
                }
                d.x=0;d.y=0;
                // Recalculating the center of each cell
                for(int h=0;h<4;h++)
                {
                    if(h==3)
                        t=0;
                    else
                        t=h+1;
                    d.x=(add(add(b[i][j][l].position[h],b[i][j][l].position[t]),b[i][j][l].center).x)*b[i][j][l].volumet[h]+d.x;
                    d.y=(add(add(b[i][j][l].position[h],b[i][j][l].position[t]),b[i][j][l].center).y)*b[i][j][l].volumet[h]+d.y;
                }
                b[i][j][l].center.x=d.x/(3*b[i][j][l].volume);
                b[i][j][l].center.y=d.y/(3*b[i][j][l].volume);
                // Calculation of area vectors bounding each cell, centre of each side of each cell, unit normal vectors for each
                // area and unit vectors for each side of the cell
                for(int h=0;h<4;h++)
                {
                    if(h==3)
                        t=0;
                    else
                        t=h+1;
                    b[i][j][l].area[h].x=b[i][j][l].position[t].y-b[i][j][l].position[h].y;
                    b[i][j][l].area[h].y=b[i][j][l].position[h].x-b[i][j][l].position[t].x;
                    b[i][j][l].side[h]=subtract(b[i][j][l].position[t],b[i][j][l].position[h]);
                    d=add(b[i][j][l].position[h],b[i][j][l].position[t]);
                    b[i][j][l].ifc[h].x=d.x/2;
                    b[i][j][l].ifc[h].y=d.y/2;
                    b[i][j][l].nn[h].x=b[i][j][l].area[h].x/mod1(b[i][j][l].area[h]);
                    b[i][j][l].nn[h].y=b[i][j][l].area[h].y/mod1(b[i][j][l].area[h]);
                    b[i][j][l].pp[h].x=b[i][j][l].side[h].x/mod1(b[i][j][l].side[h]);
                    b[i][j][l].pp[h].y=b[i][j][l].side[h].y/mod1(b[i][j][l].side[h]);
                }
                d=subtract(b[i][j][l].ifc[2],b[i][j][l].ifc[0]);
                b[i][j][l].widy=sqrt(pow(d.x,2)+pow(d.y,2));
                d=subtract(b[i][j][l].ifc[1],b[i][j][l].ifc[3]);
                b[i][j][l].widx=sqrt(pow(d.x,2)+pow(d.y,2));
            }
        }
    }
}
double roeavg(double fl,double fr,int i,int j,int k,int ll,int tt,cell***b)
{
    double mu;
    mu=(sqrt(b[i][j][ll].rhol[k])*fl+sqrt(b[i][j][tt].rhor[k])*fr)/(sqrt(b[i][j][ll].rhol[k])+sqrt(b[i][j][tt].rhor[k]));
    return mu;

}

// This function converts the velocities in cartesian coordinates into local coordinates for each interface or each bounding area of the cell
void localtocart(cell***b,cell****bound1,cell****bound2,int nx,int ny)
{
    int j,k,g;
    for(int l=0;l<2;l++)
    {
        #pragma omp parallel for private(j,k,g)
        for(int i=0;i<(nx-1);i++)
        {
            for(j=0;j<(ny-1);j++)
            {
                for(k=0;k<4;k++)
                {
                    b[i][j][l].vlocal[k].y=dot(b[i][j][l].vcart,b[i][j][l].nn[k]);
                    b[i][j][l].vlocal[k].x=dot(b[i][j][l].vcart,b[i][j][l].pp[k]);
                    for(g=0;g<4;g++)
                    {
                        bound1[i][j][l][k].vlocal[g].y=dot(bound1[i][j][l][k].vcart,bound1[i][j][l][k].nn[g]);
                        bound1[i][j][l][k].vlocal[g].x=dot(bound1[i][j][l][k].vcart,bound1[i][j][l][k].pp[g]);
                        bound2[i][j][l][k].vlocal[g].y=dot(bound2[i][j][l][k].vcart,bound2[i][j][l][k].nn[g]);
                        bound2[i][j][l][k].vlocal[g].x=dot(bound2[i][j][l][k].vcart,bound2[i][j][l][k].pp[g]);
                    }
                   // if(i>96)
                   // cout<<b[i][j][l].vlocal[k].y<<"\t"<<k<<"\n";
                }
            }
        }
    }
}
// This function calculates interfacial pressure between the two phases inside each cell
void intpressure(cell***b,double Cpstar,int nx,int ny)
{
    double deltapstar[nx-1][ny-1];
    int j;
    #pragma omp parallel for private(j)
    for(int i=0;i<(nx-1);i++)
    {
        for(j=0;j<(ny-1);j++)
        {
            deltapstar[i][j]=Cpstar*b[i][j][1].alpha*b[i][j][0].rho*pow(mod1(subtract(b[i][j][1].vcart,b[i][j][0].vcart)),2);
            deltapstar[i][j]=min(deltapstar[i][j],0.01*b[i][j][0].p);
            b[i][j][1].pint=b[i][j][1].p-deltapstar[i][j];
            b[i][j][0].pint=b[i][j][0].p-deltapstar[i][j];
        }
    }
}
// This function stores the value of conservative variables at the nth time step in an array W
void variable_calculator(cell***b,int nx,int ny)
{
    int j;
    for(int l=0;l<2;l++)
    {
        #pragma omp parallel for private(j)
        for(int i=0;i<(nx-1);i++)
        {
            for(j=0;j<(ny-1);j++)
            {
                b[i][j][l].W[0]=b[i][j][l].U[0]=b[i][j][l].alpha*b[i][j][l].rho;
                b[i][j][l].W[1]=b[i][j][l].U[1]=b[i][j][l].alpha*b[i][j][l].rho*b[i][j][l].vcart.x;
                b[i][j][l].W[2]=b[i][j][l].U[2]=b[i][j][l].alpha*b[i][j][l].rho*b[i][j][l].vcart.y;
                b[i][j][l].W[3]=b[i][j][l].U[3]=b[i][j][l].alpha*b[i][j][l].rho*b[i][j][l].E+b[i][j][l].pint*b[i][j][l].alpha;
                //cout<<b[0][0][0].U[0]<<"\t"<<b[0][0][0].U[1]<<"\t"<<b[0][0][0].U[2]<<"\t"<<b[0][0][0].U[3]<<"\n";
            }
        }
    }
}
// This function implements the boundary conditions in the sense that it stores all the values pertaining to the boundary cells of each cell
// in a new struct array bound which is again of cell type. When we reach the boundary of the control volume, we assign the boundary conditions to
// this structure bound while for an intermediate cell inside a control volume, we simply pass on the values of boundary cells to the
// array bound. Here the index 0 represents lower boundary for each cell, index 1 represents right boundary, 2 represents upper boundary
// while 3 represents left boundary
void boundary(cell***b,cell****bound1,cell****bound2,int nx,int ny,double eee)
{
    int j;
    for(int l=0;l<2;l++)
    {
        #pragma omp parallel for private(j)
        for(int i=0;i<(nx-1);i++)
        {
            for(j=0;j<(ny-1);j++)
            {
                // When j=0 or j=1, we are at the boundary of control volume, hence we implement boundary conditions
                if(j==0)
                {
                    bound1[i][j][l][0]=bound2[i][j][l][0]=b[i][j][l];
                    bound1[i][j][l][0].vcart.y=bound2[i][j][l][0].vcart.y=-b[i][j][l].vcart.y;
                }
                else if(j==1)
                {
                    bound2[i][j][l][0]=bound1[i][j][l][0]=b[i][j-1][l];
                    bound2[i][j][l][0].vcart.y=-b[i][j-1][l].vcart.y;
                }
                // when we are at an intermediate location, we simply pass on the information present in neighbouring cells to the
                // struct array bound. Similar procedure is followed for all the other boundaries.
                else
                {
                    bound1[i][j][l][0]=b[i][j-1][l];
                    bound2[i][j][l][0]=b[i][j-2][l];
                }
                if(i==(nx-2))
                {
                    bound1[i][j][l][1]=bound2[i][j][l][1]=b[i][j][l];
                    bound1[i][j][l][1].p=bound2[i][j][l][1].p=pow(10,5);
                    bound2[i][j][l][1].rho=bound1[i][j][l][1].rho=bound1[i][j][l][1].gamma*(bound1[i][j][l][1].p+bound1[i][j][l][1].pinf)/(bound1[i][j][l][1].Cp*bound1[i][j][l][1].T*(bound1[i][j][l][1].gamma-1));
                    bound2[i][j][l][1].e=bound1[i][j][l][1].e=(bound1[i][j][l][1].Cp*bound1[i][j][l][1].T/bound1[i][j][l][1].gamma)+bound1[i][j][l][1].pinf/bound1[i][j][l][1].rho;
                }
                else if(i==(nx-3))
                {
                    bound2[i][j][l][1]=bound1[i][j][l][1]=b[i+1][j][l];
                    bound2[i][j][l][1].p=pow(10,5);
                    bound2[i][j][l][1].rho=bound2[i][j][l][1].gamma*(bound2[i][j][l][1].p+bound2[i][j][l][1].pinf)/(bound2[i][j][l][1].Cp*bound2[i][j][l][1].T*(bound2[i][j][l][1].gamma-1));
                    bound2[i][j][l][1].e=(bound2[i][j][l][1].Cp*bound2[i][j][l][1].T/bound2[i][j][l][1].gamma)+bound2[i][j][l][1].pinf/bound2[i][j][l][1].rho;
                }
                else
                {
                    bound1[i][j][l][1]=b[i+1][j][l];
                    bound2[i][j][l][1]=b[i+2][j][l];
                }
                if(j==(ny-2))
                {
                    bound1[i][j][l][2]=bound2[i][j][l][2]=b[i][j][l];
                }
                else if(j==(ny-3))
                {
                    bound1[i][j][l][2]=bound2[i][j][l][2]=b[i][j+1][l];
                }
                else
                {
                    bound1[i][j][l][2]=b[i][j+1][l];
                    bound2[i][j][l][2]=b[i][j+2][l];
                }
                if(i==0)
                {
                    bound1[i][j][l][3]=bound2[i][j][l][3]=b[i][j][l];
                    //bound1[i][j][l][3].p=bound2[i][j][l][3].p=b[i][j][l].p;
                    bound1[i][j][l][3].vcart.x=bound2[i][j][l][3].vcart.x=661.81;
                    bound1[i][j][l][3].vcart.y=bound2[i][j][l][3].vcart.y=0;
                    bound1[i][j][l][3].T=bound2[i][j][l][3].T=595.13;
                    bound2[i][j][l][3].rho=bound1[i][j][l][3].rho=bound1[i][j][l][3].gamma*(bound1[i][j][l][3].p+bound1[i][j][l][3].pinf)/(bound1[i][j][l][3].Cp*bound1[i][j][l][3].T*(bound1[i][j][l][3].gamma-1));
                    bound2[i][j][l][3].e=bound1[i][j][l][3].e=(bound1[i][j][l][3].Cp*bound1[i][j][l][3].T/bound1[i][j][l][3].gamma)+bound1[i][j][l][3].pinf/bound1[i][j][l][3].rho;
                    if(l==0)
                        bound1[i][j][l][3].alpha=bound2[i][j][l][3].alpha=eee;
                    else
                        bound1[i][j][l][3].alpha=bound2[i][j][l][3].alpha=1-eee;
                }
                else if(i==1)
                {
                    bound2[i][j][l][3]=bound1[i][j][l][3]=b[i-1][j][l];
                    //bound2[i][j][l][3].p=b[i-1][j][l].p;
                    bound2[i][j][l][3].vcart.x=661.81;
                    bound2[i][j][l][3].vcart.y=0;
                    bound2[i][j][l][3].T=595.13;
                    bound2[i][j][l][3].rho=bound2[i][j][l][3].gamma*(bound2[i][j][l][3].p+bound2[i][j][l][3].pinf)/(bound2[i][j][l][3].Cp*bound2[i][j][l][3].T*(bound2[i][j][l][3].gamma-1));
                    bound2[i][j][l][3].e=(bound2[i][j][l][3].Cp*bound2[i][j][l][3].T/bound2[i][j][l][3].gamma)+bound2[i][j][l][3].pinf/bound2[i][j][l][3].rho;
                    if(l==0)
                        bound2[i][j][l][3].alpha=eee;
                    else
                        bound2[i][j][l][3].alpha=1-eee;
                }
                else
                {
                    bound1[i][j][l][3]=b[i-1][j][l];
                    bound2[i][j][l][3]=b[i-2][j][l];
                }
            }
        }
    }
}
// This function is used to reconstruct the values of variables at the cell boundaries. We use the same indices for identifying cell
// boundaries as mentioned earlier.
// We reconstruct the variables using a one dimensional approach for each cell interface. bound1 and bound2 hold the values of variables
// for first two neighbouring cells in the corresponding direction.
// We only reconstruct velocity, internal energy, volume fraction and density. Rest are found using equation of state and other equations
void reconstruction(cell***b,cell****bound1,cell****bound2,int gg,int nx,int ny)
{
    double ee,alphal,alphar;
    ee=pow(10,-12);
    int i,k,ct;
    double***l0=new double**[ny-1];
    double***r0=new double**[ny-1];
    double***l1=new double**[ny-1];
    double***r1=new double**[ny-1];
    double**deltaln=new double*[ny-1];
    double**deltarn=new double*[ny-1];
    double**deltarp=new double*[ny-1];
    double**deltalp=new double*[ny-1];
    double**sl=new double*[ny-1];
    double**sr=new double*[ny-1];
    double**h=new double*[ny-1];
    for(int j=0;j<(ny-1);j++)
    {
        l0[j]=new double*[4];
        l1[j]=new double*[4];
        r0[j]=new double*[4];
        r1[j]=new double*[4];
        deltaln[j]=new double[5];
        deltarn[j]=new double[5];
        deltarp[j]=new double[5];
        deltalp[j]=new double[5];
        sl[j]=new double[5];
        sr[j]=new double[5];
        h[j]=new double[4];
        for(i=0;i<4;i++)
        {
            l0[j][i]=new double[5];
            l1[j][i]=new double[5];
            r0[j][i]=new double[5];
            r1[j][i]=new double[5];
        }
    }
    int j;
    // Reconstruction procedure starts here. For details of procedure and equations, refer Eilmer.
    for(int l=0;l<2;l++)
    {
        #pragma omp parallel for private(j,k,ct,alphal,alphar)
        for(int i=0;i<(nx-1);i++)
        {
            for(j=0;j<(ny-1);j++)
            {
                // Reconstruction takes place for each interface
                for(k=0;k<4;k++)
                {
                    // l0 stores values of quantities to be reconstructed of the cell in consideration
                    l0[j][k][0]=b[i][j][l].rho;
                    l0[j][k][1]=b[i][j][l].vlocal[k].x;
                    l0[j][k][2]=b[i][j][l].vlocal[k].y;
                    l0[j][k][3]=b[i][j][l].e;
                    l0[j][k][4]=b[i][j][l].alpha;
                    // if(k==0), we reconstruct the values for the bottom interface
                    if(k==0)
                    {
                         l1[j][k][0]=bound1[i][j][l][2].rho;
                         l1[j][k][1]=bound1[i][j][l][2].vlocal[k].x;
                         l1[j][k][2]=bound1[i][j][l][2].vlocal[k].y;
                         l1[j][k][3]=bound1[i][j][l][2].e;
                         l1[j][k][4]=bound1[i][j][l][2].alpha;
                         r0[j][k][0]=bound1[i][j][l][0].rho;
                         r0[j][k][1]=bound1[i][j][l][0].vlocal[k].x;
                         r0[j][k][2]=bound1[i][j][l][0].vlocal[k].y;
                         r0[j][k][3]=bound1[i][j][l][0].e;
                         r0[j][k][4]=bound1[i][j][l][0].alpha;
                         r1[j][k][0]=bound2[i][j][l][0].rho;
                         r1[j][k][1]=bound2[i][j][l][0].vlocal[k].x;
                         r1[j][k][2]=bound2[i][j][l][0].vlocal[k].y;
                         r1[j][k][3]=bound2[i][j][l][0].e;
                         r1[j][k][4]=bound2[i][j][l][0].alpha;
                         h[j][0]=bound1[i][j][l][2].widy;
                         h[j][1]=b[i][j][l].widy;
                         h[j][2]=bound1[i][j][l][0].widy;
                         h[j][3]=bound2[i][j][l][0].widy;
                    }
                    // if(k==1), we reconstruct the values for the right interface
                    else if(k==1)
                    {
                        l1[j][k][0]=bound1[i][j][l][3].rho;
                        l1[j][k][1]=bound1[i][j][l][3].vlocal[k].x;
                        l1[j][k][2]=bound1[i][j][l][3].vlocal[k].y;
                        l1[j][k][3]=bound1[i][j][l][3].e;
                        l1[j][k][4]=bound1[i][j][l][3].alpha;
                        r0[j][k][0]=bound1[i][j][l][1].rho;
                        r0[j][k][1]=bound1[i][j][l][1].vlocal[k].x;
                        r0[j][k][2]=bound1[i][j][l][1].vlocal[k].y;
                        r0[j][k][3]=bound1[i][j][l][1].e;
                        r0[j][k][4]=bound1[i][j][l][1].alpha;
                        r1[j][k][0]=bound2[i][j][l][1].rho;
                        r1[j][k][1]=bound2[i][j][l][1].vlocal[k].x;
                        r1[j][k][2]=bound2[i][j][l][1].vlocal[k].y;
                        r1[j][k][3]=bound2[i][j][l][1].e;
                        r1[j][k][4]=bound2[i][j][l][1].alpha;
                        h[j][0]=bound1[i][j][l][3].widx;
                        h[j][1]=b[i][j][l].widx;
                        h[j][2]=bound1[i][j][l][1].widx;
                        h[j][3]=bound2[i][j][l][1].widx;
                    }
                    // if(k==2), we reconstruct the values for the upper interface
                    else if(k==2)
                    {
                        l1[j][k][0]=bound1[i][j][l][0].rho;
                        l1[j][k][1]=bound1[i][j][l][0].vlocal[k].x;
                        l1[j][k][2]=bound1[i][j][l][0].vlocal[k].y;
                        l1[j][k][3]=bound1[i][j][l][0].e;
                        l1[j][k][4]=bound1[i][j][l][0].alpha;
                        r0[j][k][0]=bound1[i][j][l][2].rho;
                        r0[j][k][1]=bound1[i][j][l][2].vlocal[k].x;
                        r0[j][k][2]=bound1[i][j][l][2].vlocal[k].y;
                        r0[j][k][3]=bound1[i][j][l][2].e;
                        r0[j][k][4]=bound1[i][j][l][2].alpha;
                        r1[j][k][0]=bound2[i][j][l][2].rho;
                        r1[j][k][1]=bound2[i][j][l][2].vlocal[k].x;
                        r1[j][k][2]=bound2[i][j][l][2].vlocal[k].y;
                        r1[j][k][3]=bound2[i][j][l][2].e;
                        r1[j][k][4]=bound2[i][j][l][2].alpha;
                        h[j][0]=bound1[i][j][l][0].widy;
                        h[j][1]=b[i][j][l].widy;
                        h[j][2]=bound1[i][j][l][2].widy;
                        h[j][3]=bound2[i][j][l][2].widy;
                    }
                    // if(k==3), we reconstruct the values for the left interface
                    else if(k==3)
                    {
                        l1[j][k][0]=bound1[i][j][l][1].rho;
                        l1[j][k][1]=bound1[i][j][l][1].vlocal[k].x;
                        l1[j][k][2]=bound1[i][j][l][1].vlocal[k].y;
                        l1[j][k][3]=bound1[i][j][l][1].e;
                        l1[j][k][4]=bound1[i][j][l][1].alpha;
                        r0[j][k][0]=bound1[i][j][l][3].rho;
                        r0[j][k][1]=bound1[i][j][l][3].vlocal[k].x;
                        r0[j][k][2]=bound1[i][j][l][3].vlocal[k].y;
                        r0[j][k][3]=bound1[i][j][l][3].e;
                        r0[j][k][4]=bound1[i][j][l][3].alpha;
                        r1[j][k][0]=bound2[i][j][l][3].rho;
                        r1[j][k][1]=bound2[i][j][l][3].vlocal[k].x;
                        r1[j][k][2]=bound2[i][j][l][3].vlocal[k].y;
                        r1[j][k][3]=bound2[i][j][l][3].e;
                        r1[j][k][4]=bound1[i][j][l][3].alpha;
                        h[j][0]=bound1[i][j][l][1].widx;
                        h[j][1]=b[i][j][l].widx;
                        h[j][2]=bound1[i][j][l][3].widx;
                        h[j][3]=bound2[i][j][l][3].widx;
                    }
                    // Refer eilmer for the description of equations and variables used below
                    for(ct=0;ct<5;ct++)
                    {
                        deltaln[j][ct]=2*(l0[j][k][ct]-l1[j][k][ct])/(h[j][1]+h[j][0]);
                        deltalp[j][ct]=deltarn[j][ct]=2*(r0[j][k][ct]-l0[j][k][ct])/(h[j][2]+h[j][1]);
                        deltarp[j][ct]=2*(r1[j][k][ct]-r0[j][k][ct])/(h[j][2]+h[j][3]);
                        sl[j][ct]=(deltaln[j][ct]*deltalp[j][ct]+fabs(deltaln[j][ct]*deltalp[j][ct]))/(pow(deltaln[j][ct],2)+pow(deltalp[j][ct],2)+ee);
                        sr[j][ct]=(deltarn[j][ct]*deltarp[j][ct]+fabs(deltarn[j][ct]*deltarp[j][ct]))/(pow(deltarn[j][ct],2)+pow(deltarp[j][ct],2)+ee);
                    }
                    alphal=h[j][1]/(2*(h[j][0]+2*h[j][1]+h[j][2]));
                    alphar=h[j][2]/(2*(h[j][1]+2*h[j][2]+h[j][3]));

                    // Reconstruction of density, velocity, internal energy, volume fraction
                    b[i][j][l].rhol[k]=l0[j][k][0]+alphal*(deltalp[j][0]*(2*h[j][1]+h[j][0])+deltaln[j][0]*h[j][2])*sl[j][0];
                    b[i][j][l].rhor[k]=r0[j][k][0]-alphar*(deltarn[j][0]*(2*h[j][2]+h[j][3])+deltarp[j][0]*h[j][1])*sr[j][0];
                    b[i][j][l].vllocal[k].x=l0[j][k][1]+alphal*(deltalp[j][1]*(2*h[j][1]+h[j][0])+deltaln[j][1]*h[j][2])*sl[j][1];
                    b[i][j][l].vrlocal[k].x=r0[j][k][1]-alphar*(deltarn[j][1]*(2*h[j][2]+h[j][3])+deltarp[j][1]*h[j][1])*sr[j][1];
                    b[i][j][l].vllocal[k].y=l0[j][k][2]+alphal*(deltalp[j][2]*(2*h[j][1]+h[j][0])+deltaln[j][2]*h[j][2])*sl[j][2];
                    b[i][j][l].vrlocal[k].y=r0[j][k][2]-alphar*(deltarn[j][2]*(2*h[j][2]+h[j][3])+deltarp[j][2]*h[j][1])*sr[j][2];
                    b[i][j][l].el[k]=l0[j][k][3]+alphal*(deltalp[j][3]*(2*h[j][1]+h[j][0])+deltaln[j][3]*h[j][2])*sl[j][3];
                    b[i][j][l].er[k]=r0[j][k][3]-alphar*(deltarn[j][3]*(2*h[j][2]+h[j][3])+deltarp[j][3]*h[j][1])*sr[j][3];
                    b[i][j][l].alphal[k]=l0[j][k][4]+alphal*(deltalp[j][4]*(2*h[j][1]+h[j][0])+deltaln[j][4]*h[j][2])*sl[j][4];
                    b[i][j][l].alphar[k]=r0[j][k][4]-alphar*(deltarn[j][4]*(2*h[j][2]+h[j][3])+deltarp[j][4]*h[j][1])*sr[j][4];

                    // We use the min-max limiter to avoid large variations in reconstructed quantities from the cell values
                    if(b[i][j][l].rhol[k]<min(l0[j][k][0],r0[j][k][0]))
                        b[i][j][l].rhol[k]=min(l0[j][k][0],r0[j][k][0]);
                    else if(b[i][j][l].rhol[k]>max(l0[j][k][0],r0[j][k][0]))
                        b[i][j][l].rhol[k]=max(l0[j][k][0],r0[j][k][0]);
                    if(b[i][j][l].rhor[k]<min(l0[j][k][0],r0[j][k][0]))
                        b[i][j][l].rhor[k]=min(l0[j][k][0],r0[j][k][0]);
                    else if(b[i][j][l].rhor[k]>max(l0[j][k][0],r0[j][k][0]))
                        b[i][j][l].rhor[k]=max(l0[j][k][0],r0[j][k][0]);

                    if(b[i][j][l].vllocal[k].x<min(l0[j][k][1],r0[j][k][1]))
                        b[i][j][l].vllocal[k].x=min(l0[j][k][1],r0[j][k][1]);
                    else if(b[i][j][l].vllocal[k].x>max(l0[j][k][1],r0[j][k][1]))
                        b[i][j][l].vllocal[k].x=max(l0[j][k][1],r0[j][k][1]);
                    if(b[i][j][l].vrlocal[k].x<min(l0[j][k][1],r0[j][k][1]))
                        b[i][j][l].vrlocal[k].x=min(l0[j][k][1],r0[j][k][1]);
                    else if(b[i][j][l].vrlocal[k].x>max(l0[j][k][1],r0[j][k][1]))
                        b[i][j][l].vrlocal[k].x=max(l0[j][k][1],r0[j][k][1]);

                    if(b[i][j][l].vllocal[k].y<min(l0[j][k][2],r0[j][k][2]))
                        b[i][j][l].vllocal[k].y=min(l0[j][k][2],r0[j][k][2]);
                    else if(b[i][j][l].vllocal[k].y>max(l0[j][k][2],r0[j][k][2]))
                        b[i][j][l].vllocal[k].y=max(l0[j][k][2],r0[j][k][2]);
                    if(b[i][j][l].vrlocal[k].y<min(l0[j][k][2],r0[j][k][2]))
                        b[i][j][l].vrlocal[k].y=min(l0[j][k][2],r0[j][k][2]);
                    else if(b[i][j][l].vrlocal[k].y>max(l0[j][k][2],r0[j][k][2]))
                        b[i][j][l].vrlocal[k].y=max(l0[j][k][2],r0[j][k][2]);

                    if(b[i][j][l].el[k]<min(l0[j][k][3],r0[j][k][3]))
                        b[i][j][l].el[k]=min(l0[j][k][3],r0[j][k][3]);
                    else if(b[i][j][l].el[k]>max(l0[j][k][3],r0[j][k][3]))
                        b[i][j][l].el[k]=max(l0[j][k][3],r0[j][k][3]);
                    if(b[i][j][l].er[k]<min(l0[j][k][3],r0[j][k][3]))
                        b[i][j][l].er[k]=min(l0[j][k][3],r0[j][k][3]);
                    else if(b[i][j][l].er[k]>max(l0[j][k][3],r0[j][k][3]))
                        b[i][j][l].er[k]=max(l0[j][k][3],r0[j][k][3]);

                    if(b[i][j][l].alphal[k]<min(l0[j][k][4],r0[j][k][4]))
                        b[i][j][l].alphal[k]=min(l0[j][k][4],r0[j][k][4]);
                    else if(b[i][j][l].alphal[k]>max(l0[j][k][4],r0[j][k][4]))
                        b[i][j][l].alphal[k]=max(l0[j][k][4],r0[j][k][4]);
                    if(b[i][j][l].alphar[k]<min(l0[j][k][4],r0[j][k][4]))
                        b[i][j][l].alphar[k]=min(l0[j][k][4],r0[j][k][4]);
                    else if(b[i][j][l].alphar[k]>max(l0[j][k][4],r0[j][k][4]))
                        b[i][j][l].alphar[k]=max(l0[j][k][4],r0[j][k][4]);

                    b[i][j][l].pl[k]=b[i][j][l].rhol[k]*b[i][j][l].el[k]*(b[i][j][l].gamma-1)-b[i][j][l].gamma*b[i][j][l].pinf;
                    b[i][j][l].pr[k]=b[i][j][l].rhor[k]*b[i][j][l].er[k]*(b[i][j][l].gamma-1)-b[i][j][l].gamma*b[i][j][l].pinf;

                    if(b[i][j][l].pl[k]<min(b[i][j][l].p,bound1[i][j][l][k].p))
                       b[i][j][l].pl[k]=min(b[i][j][l].p,bound1[i][j][l][k].p);
                    else if(b[i][j][l].pl[k]>max(b[i][j][l].p,bound1[i][j][l][k].p))
                       b[i][j][l].pl[k]=max(b[i][j][l].p,bound1[i][j][l][k].p);
                    if(b[i][j][l].pr[k]<min(b[i][j][l].p,bound1[i][j][l][k].p))
                       b[i][j][l].pr[k]=min(b[i][j][l].p,bound1[i][j][l][k].p);
                    else if(b[i][j][l].pr[k]>max(b[i][j][l].p,bound1[i][j][l][k].p))
                       b[i][j][l].pr[k]=max(b[i][j][l].p,bound1[i][j][l][k].p);

                    b[i][j][l].Tl[k]=(b[i][j][l].pl[k]+b[i][j][l].pinf)*b[i][j][l].gamma/((b[i][j][l].gamma-1)*b[i][j][l].rhol[k]*b[i][j][l].Cp);
                    b[i][j][l].Tr[k]=(b[i][j][l].pr[k]+b[i][j][l].pinf)*b[i][j][l].gamma/((b[i][j][l].gamma-1)*b[i][j][l].rhor[k]*b[i][j][l].Cp);

                    if(b[i][j][l].Tl[k]<min(b[i][j][l].T,bound1[i][j][l][k].T))
                       b[i][j][l].Tl[k]=min(b[i][j][l].T,bound1[i][j][l][k].T);
                    else if(b[i][j][l].Tl[k]>max(b[i][j][l].T,bound1[i][j][l][k].T))
                       b[i][j][l].Tl[k]=max(b[i][j][l].T,bound1[i][j][l][k].T);
                    if(b[i][j][l].Tr[k]<min(b[i][j][l].T,bound1[i][j][l][k].T))
                       b[i][j][l].Tr[k]=min(b[i][j][l].T,bound1[i][j][l][k].T);
                    else if(b[i][j][l].Tr[k]>max(b[i][j][l].T,bound1[i][j][l][k].T))
                       b[i][j][l].Tr[k]=max(b[i][j][l].T,bound1[i][j][l][k].T);

                    b[i][j][l].al[k]=sqrt(b[i][j][l].gamma*(b[i][j][l].pl[k]+b[i][j][l].pinf)/b[i][j][l].rhol[k]);
                    b[i][j][l].ar[k]=sqrt(b[i][j][l].gamma*(b[i][j][l].pr[k]+b[i][j][l].pinf)/b[i][j][l].rhor[k]);

                    b[i][j][l].Hl[k]=b[i][j][l].el[k]+pow(mod1(b[i][j][l].vllocal[k]),2)/2+b[i][j][l].pl[k]/b[i][j][l].rhol[k];
                    b[i][j][l].Hr[k]=b[i][j][l].er[k]+pow(mod1(b[i][j][l].vrlocal[k]),2)/2+b[i][j][l].pr[k]/b[i][j][l].rhor[k];
                    //if(l==1)
                     //cout<<b[0][0][l].vllocal[k].y<<"\t"<<b[0][0][l].vlocal[k].y<<"\t"<<b[0][0][l].vrlocal[k].y<<"\t"<<k<<"\t"<<l<<"\n";

                    // deltaeff stores the difference between the left and the right gas volume fractions for each interface
                    // In other words, deltaeff represents the extent of phase transition
                    b[i][j][l].deltaeff[k]=fabs(b[i][j][0].alphal[k]-b[i][j][0].alphar[k]);
                }
            }
        }
    }
    for(int j=0;j<(ny-1);j++)
    {
        for(int i=0;i<4;i++)
        {
            delete[] l0[j][i];
            delete[] r0[j][i];
            delete[] l1[j][i];
            delete[] r1[j][i];
        }
        delete[] l0[j];
        delete[] r0[j];
        delete[] l1[j];
        delete[] r1[j];
        delete[] deltaln[j];
        delete[] deltalp[j];
        delete[] deltarn[j];
        delete[] deltarp[j];
        delete[] sl[j];
        delete[] sr[j];
        delete[] h[j];
    }
}

// Function calculates speed of sound at the boundary separating two cells
void intsoundspeed(cell***b,int nx,int ny,double eee)
{
    int j,k;
    for(int l=0;l<2;l++)
    {
        #pragma omp parallel for private(j,k)
        for(int i=0;i<(nx-1);i++)
        {
            for(j=0;j<(ny-1);j++)
            {
                for(k=0;k<4;k++)
                {
                    b[i][j][l].ic_ak[k]=(b[i][j][l].al[k]+b[i][j][l].ar[k])/2;

                }
            }
        }
    }
    for(int l=0;l<2;l++)
    {
        #pragma omp parallel for private(j,k)
        for(int i=0;i<(nx-1);i++)
        {
            for(j=0;j<(ny-1);j++)
            {
                for(k=0;k<4;k++)
                {
                   /* if(b[i][j][l].deltaeff[k]<5*eee)
                    {
                        b[i][j][l].ic_a[k]=(b[i][j][0].ic_ak[k]+b[i][j][1].ic_ak[k])/2;
                    }
                    else*/
                    {
                        b[i][j][l].ic_a[k]=b[i][j][l].ic_ak[k];
                    }
                    //cout<<b[i][j][l].ic_a[k]<<"\n";
                }
            }
        }
    }
}

// This function calculates the flux terms. Here we first calculate the gas-gas and liquid-liquid flux at each boundary using AUSM+-up scheme
// Then using HLLC scheme, we calculate the gas-liquid or liquid-gas flux at each boundary.
// While calculating flux,we first calculate them in the local frame of reference for each interface. Then each flux term is converted into
// global coordinates and added up to find the total flux term
// Another thing to be noted is that it is only when the phase transition between neighbouring cells is greater than a specified
// tolerance value that we use both HLLC and AUSM+-up scheme. In case, phase transition is less than a specifed tolerance value, we use only
// the AUSM+-up scheme to calculate the total flux term
void interfacial(cell***b,int nx,int ny, double Kp, double Ku, double eee)
{
    double ic_rho,Mbar,Mp,Pu,M1positive,M2positive,M4positive,M1negative,M2negative,M4negative; // Refer AUSM+-up scheme
    double P5positive,P5negative,ML,MR,mdot,ic_M,ic_p; // Refer AUSM+-up scheme
    int i,k,jj;
    // ic_fluxc stores convective flux term, ic_fluxp stores pressure flux term in local frame of reference of each interface
    // ic_flux stores the sum of pressure and convective flux terms in local frame of reference
    // ncf stores non conservative flux term in local frame of reference
    // ic_fluxcart stores total flux in global frame of reference
    double**ic_fluxc=new double*[ny-1];
    double**ic_fluxp=new double*[ny-1];
    double**ic_flux=new double*[ny-1];
    double**ncf=new double*[ny-1];
    double**ic_fluxn=new double*[ny-1];
    double**ic_fluxcart=new double*[ny-1];
    for(i=0;i<(ny-1);i++)
    {
        ic_fluxc[i]=new double[4];
        ic_fluxp[i]=new double[4];
        ic_flux[i]=new double[4];
        ncf[i]=new double[4];
        ic_fluxn[i]=new double[4];
        ic_fluxcart[i]=new double[4];
    }
    int j;
    // Calculation of gas-gas or liquid-liquid flux using AUSM+-up scheme
    for(int l=0;l<2;l++)
    {
        #pragma omp parallel for private(j,k,jj,ic_rho,Mbar,Mp,Pu,M1positive,M2positive,M4positive,M1negative,M2negative,M4negative,P5positive,P5negative,ML,MR,mdot,ic_M,ic_p)
        for(int i=0;i<(nx-1);i++)
        {
            for(j=0;j<(ny-1);j++)
            {
                // for each cell, we assign initial values to the total flux r as 0;
                b[i][j][l].r[0]=0;b[i][j][l].r[1]=0;b[i][j][l].r[2]=0;b[i][j][l].r[3]=0;
                for(k=0;k<4;k++)
                {
                    // Refer AUSM+-up scheme. usual process of evaluation of variables.
                    ML=b[i][j][l].vllocal[k].y/b[i][j][l].ic_a[k];
                    MR=b[i][j][l].vrlocal[k].y/b[i][j][l].ic_a[k];
                    ic_rho=(b[i][j][l].rhol[k]+b[i][j][l].rhor[k])/2;
                    Mbar=sqrt((pow(ML,2)+pow(MR,2))/2);
                    Mp=-Kp*max(1-pow(Mbar,2),0.0)*(b[i][j][l].pr[k]-b[i][j][l].pl[k])/(ic_rho*pow(b[i][j][l].ic_a[k],2));
                    M1positive=0.5*(ML+fabs(ML));
                    M2positive=0.25*pow((ML+1),2);
                    M1negative=0.5*(MR-fabs(MR));
                    M2negative=-0.25*pow(MR-1,2);
                    if(fabs(ML)>1)
                    {
                        M4positive=M1positive;
                        P5positive=M1positive/ML;
                    }
                    else
                    {
                        M4positive=M2positive*(1-2*M2negative);
                        P5positive=M2positive*(2-ML-3*ML*M2negative);
                    }
                    if(fabs(MR)>1)
                    {
                        M4negative=M1negative;
                        P5negative=M1negative/MR;
                    }
                    else
                    {
                        M4negative=M2negative*(1+2*M2positive);
                        P5negative=M2negative*(-2-MR+3*MR*M2positive);
                    }
                    Pu=-Ku*P5positive*P5negative*ic_rho*b[i][j][l].ic_a[k]*(b[i][j][l].vrlocal[k].y-b[i][j][l].vllocal[k].y);
                    ic_M=M4positive+M4negative+Mp;
                    if(ic_M>0)
                    {
                        mdot=ic_M*b[i][j][l].ic_a[k]*b[i][j][l].rhol[k];
                    }
                    else
                    {
                        mdot=ic_M*b[i][j][l].ic_a[k]*b[i][j][l].rhor[k];
                    }
                    //cout<<mdot<<"\n";
                    ic_p=P5positive*b[i][j][l].pl[k]+P5negative*b[i][j][l].pr[k]+Pu;
                    // Evaluation of convective flux term
                    if(mdot>0)
                    {
                        ic_fluxc[j][0]=mdot*b[i][j][l].alphal[k];
                        ic_fluxc[j][1]=mdot*b[i][j][l].alphal[k]*b[i][j][l].vllocal[k].y;
                        ic_fluxc[j][2]=mdot*b[i][j][l].alphal[k]*b[i][j][l].vllocal[k].x;
                        ic_fluxc[j][3]=mdot*b[i][j][l].alphal[k]*b[i][j][l].Hl[k];
                    }
                    else
                    {
                        ic_fluxc[j][0]=mdot*b[i][j][l].alphar[k];
                        ic_fluxc[j][1]=mdot*b[i][j][l].alphar[k]*b[i][j][l].vrlocal[k].y;
                        ic_fluxc[j][2]=mdot*b[i][j][l].alphar[k]*b[i][j][l].vrlocal[k].x;
                        ic_fluxc[j][3]=mdot*b[i][j][l].alphar[k]*b[i][j][l].Hr[k];
                    }
                    ic_fluxp[j][0]=0;
                    // evaluation of pressure flux.
                    if(b[i][j][0].deltaeff[k]>eee)
                    {
                        ic_fluxp[j][1]=ic_p*min(b[i][j][l].alphal[k],b[i][j][l].alphar[k]);
                    }
                    else
                    {
                        ic_fluxp[j][1]=ic_p*b[i][j][l].alphal[k];
                    }
                    ic_fluxp[j][2]=0;
                    ic_fluxp[j][3]=0;
                    //cout<<ic_fluxp[1]<<"\t"<<i<<"\t"<<j<<"\t"<<"\n";
                    // Evaluation of sum of pressure and convective flux term
                    for(jj=0;jj<4;jj++)
                    {
                        ic_flux[j][jj]=(ic_fluxp[j][jj]+ic_fluxc[j][jj])*mod1(b[i][j][l].area[k]);
                    }
                    //  Evaluation of non conservative flux term
                    ncf[j][0]=0;
                    ncf[j][1]=(b[i][j][l].alphal[k]*b[i][j][l].pint)*mod1(b[i][j][l].area[k]);
                    ncf[j][2]=0;
                    ncf[j][3]=0;
                    ic_fluxn[j][0]=-ic_flux[j][0]+ncf[j][0];
                    ic_fluxn[j][1]=-ic_flux[j][1]+ncf[j][1];
                    ic_fluxn[j][2]=-ic_flux[j][2]+ncf[j][2];
                    ic_fluxn[j][3]=-ic_flux[j][3]+ncf[j][3];
                    // Conversion of total flux term from local frame of reference to global frame of reference
                    ic_fluxcart[j][0]=ic_fluxn[j][0];
                    ic_fluxcart[j][1]=(ic_fluxn[j][1]*b[i][j][l].nn[k].x+ic_fluxn[j][2]*b[i][j][l].pp[k].x);
                    ic_fluxcart[j][2]=(ic_fluxn[j][1]*b[i][j][l].nn[k].y+ic_fluxn[j][2]*b[i][j][l].pp[k].y);
                    ic_fluxcart[j][3]=ic_fluxn[j][3];
                    // Summing or evaluating the complete flux using AUSM+-up scheme for the entire cell. Here we add the value of
                    // ic_fluxcart for each bounding interface and calculate this term.
                    for(jj=0;jj<4;jj++)
                    {
                        b[i][j][l].r[jj]=b[i][j][l].r[jj]+ic_fluxcart[j][jj];
                    }
                   // cout<<b[i][j][l].r[0]<<"\t"<<i<<"\t"<<j<<"\t"<<"\n";
                }
            }
        }
    }
    for(int j=0;j<(ny-1);j++)
    {
        delete[] ic_fluxc[j];
        delete[] ic_fluxp[j];
        delete[] ic_flux[j];
        delete[] ic_fluxcart[j];
        delete[] ic_fluxn[j];
        delete[] ncf[j];
    }
    delete ic_fluxc;
    delete ic_fluxp;
    delete ic_flux;
    delete ic_fluxcart;
    delete ic_fluxn;
    delete ncf;

    // Evaluation of flux term for gas-liquid or liquid-gas phase transition at the boundaries using HLLC scheme. Refer to the paper by Kitamura et al for more information on equations and variables
    double vstar,pstar,rhostarl,rhostarr,bl,br,vup,gammaup,aup,pup,prhoup,fl,fr,thetal,thetar,eel,eer,hstar,alphainti,alphaintip,pint;
    // vup, aup, pup, gammaup are roe averaged variables, evaluated using the roeaverage function above.
    int ll,tt;
    // Flocal and Fcart are the fluxes for each interface in the local and global frame of reference respectively
    double**Flocal=new double*[ny-1];
    double**Fcart=new double*[ny-1];
    for(int j=0;j<(ny-1);j++)
    {
        Flocal[j]=new double[4];
        Fcart[j]=new double[4];
    }
    #pragma omp parallel for private(j,k,jj,vstar,pstar,rhostarl,rhostarr,bl,br,vup,gammaup,aup,pup,prhoup,fl,fr,thetal,thetar,eel,eer,hstar,alphainti,alphaintip,pint,ll,tt)
    for(int i=0;i<(nx-1);i++)
    {
        for(j=0;j<(ny-1);j++)
        {
            for(k=0;k<4;k++)
            {
                // We use HLLC scheme only when the phase transition at the interface is greater than a specified tolerance value
                if(b[i][j][0].deltaeff[k]<5*eee)
                {
                    continue;
                }
                else
                {
                    // We identify the dominant phase on the left and right side of the phase interface
                    if((b[i][j][0].alphal[k]+5*eee)<b[i][j][0].alphar[k])
                    {
                        ll=1;
                        tt=0;
                    }
                    else
                    {
                        ll=0;
                        tt=1;
                    }
                    // Calculation of roe averaged variables
                    fl=(b[i][j][ll].pl[k]+b[i][j][ll].pinf)/b[i][j][ll].rhol[k];
                    fr=(b[i][j][tt].pr[k]+b[i][j][tt].pinf)/b[i][j][tt].rhor[k];
                    prhoup=roeavg(fl,fr,i,j,k,ll,tt,b);
                    fl=b[i][j][ll].gamma;
                    fr=b[i][j][tt].gamma;
                    gammaup=roeavg(fl,fr,i,j,k,ll,tt,b);
                    fl=b[i][j][ll].vllocal[k].y;
                    fr=b[i][j][tt].vrlocal[k].y;
                    vup=roeavg(fl,fr,i,j,k,ll,tt,b);
                    pup=((gammaup-1)/gammaup)*(prhoup+0.5*sqrt(b[i][j][ll].rhol[k]*b[i][j][tt].rhor[k])*pow(b[i][j][ll].vllocal[k].y-b[i][j][tt].vrlocal[k].y,2)/pow(sqrt(b[i][j][ll].rhol[k])+sqrt(b[i][j][tt].rhor[k]),2));
                    aup=sqrt(gammaup*pup);
                    bl=min(b[i][j][ll].vllocal[k].y-b[i][j][ll].al[k],min((vup-aup),0.0));
                    br=max(b[i][j][tt].vrlocal[k].y+b[i][j][tt].ar[k],max(vup+aup,0.0));
                    vstar=(b[i][j][ll].rhol[k]*b[i][j][ll].vllocal[k].y*(b[i][j][ll].vllocal[k].y-bl)+b[i][j][tt].rhor[k]*b[i][j][tt].vrlocal[k].y*(br-b[i][j][tt].vrlocal[k].y)+b[i][j][ll].pl[k]-b[i][j][tt].pr[k])/(b[i][j][ll].rhol[k]*(b[i][j][ll].vllocal[k].y-bl)+b[i][j][tt].rhor[k]*(br-b[i][j][tt].vrlocal[k].y));
                    pstar=0.5*(b[i][j][ll].pl[k]+b[i][j][ll].rhol[k]*(b[i][j][ll].vllocal[k].y-bl)*(b[i][j][ll].vllocal[k].y-vstar)+b[i][j][tt].pr[k]+b[i][j][tt].rhor[k]*(br-b[i][j][tt].vrlocal[k].y)*(vstar-b[i][j][tt].vrlocal[k].y));
                    pstar=max(pstar,pow(10,-6)*min(b[i][j][ll].pl[k],b[i][j][tt].pr[k]));
                    thetal=(b[i][j][ll].gamma+1)/(b[i][j][ll].gamma-1);
                    thetar=(b[i][j][tt].gamma+1)/(b[i][j][tt].gamma-1);
                    eel=(pstar+b[i][j][ll].pinf)/(b[i][j][ll].pl[k]+b[i][j][ll].pinf);
                    eer=(pstar+b[i][j][tt].pinf)/(b[i][j][tt].pr[k]+b[i][j][tt].pinf);
                    if(pstar>b[i][j][ll].pl[k])
                    {
                        rhostarl=b[i][j][ll].rhol[k]*(thetal*eel+1)/(thetal+eel);
                    }
                    else
                    {
                        rhostarl=b[i][j][ll].rhol[k]*pow(eel,1/b[i][j][ll].gamma);
                    }
                    if(pstar>b[i][j][tt].pr[k])
                    {
                        rhostarr=b[i][j][tt].rhor[k]*(thetar*eer+1)/(thetar+eer);
                    }
                    else
                    {
                        rhostarr=b[i][j][tt].rhor[k]*pow(eer,1/b[i][j][tt].gamma);
                    }
                    // Calculation of local fluxes
                    if(vstar>0)
                    {
                        hstar=(b[i][j][ll].gamma/(b[i][j][ll].gamma-1))*(pstar+b[i][j][ll].pinf)/(rhostarl)+0.5*(pow(vstar,2)+pow(b[i][j][ll].vllocal[k].x,2));
                        Flocal[j][0]=b[i][j][0].deltaeff[k]*rhostarl*vstar;
                        Flocal[j][1]=b[i][j][0].deltaeff[k]*rhostarl*pow(vstar,2);
                        Flocal[j][2]=b[i][j][0].deltaeff[k]*rhostarl*vstar*b[i][j][ll].vllocal[k].x;
                        Flocal[j][3]=b[i][j][0].deltaeff[k]*rhostarl*vstar*hstar;
                    }
                    else
                    {
                        hstar=(b[i][j][tt].gamma/(b[i][j][tt].gamma-1))*(pstar+b[i][j][tt].pinf)/(rhostarr)+0.5*(pow(vstar,2)+pow(b[i][j][tt].vrlocal[k].x,2));
                        Flocal[j][0]=b[i][j][0].deltaeff[k]*rhostarr*vstar;
                        Flocal[j][1]=b[i][j][0].deltaeff[k]*rhostarr*pow(vstar,2);
                        Flocal[j][2]=b[i][j][0].deltaeff[k]*rhostarr*vstar*b[i][j][tt].vrlocal[k].x;
                        Flocal[j][3]=b[i][j][0].deltaeff[k]*rhostarr*vstar*hstar;
                    }
                    // Calculation of fluxes in global coordinate system
                    Fcart[j][0]=Flocal[j][0]*mod1(b[i][j][0].area[k]);
                    Fcart[j][1]=(Flocal[j][1]*b[i][j][0].nn[k].x+Flocal[j][2]*b[i][j][0].pp[k].x)*mod1(b[i][j][0].area[k]);
                    Fcart[j][2]=(Flocal[j][1]*b[i][j][0].nn[k].y+Flocal[j][2]*b[i][j][0].pp[k].y)*mod1(b[i][j][0].area[k]);
                    Fcart[j][3]=Flocal[j][3]*mod1(b[i][j][0].area[k]);
                    alphainti=0.5*(b[i][j][ll].alphal[k]+b[i][j][ll].alphar[k]);
                    alphaintip=1-alphainti;
                    pint=pstar-2*alphainti*alphaintip*b[i][j][ll].rhol[k]*b[i][j][tt].rhor[k]*pow(b[i][j][ll].vllocal[k].y-b[i][j][tt].vrlocal[k].y,2)/(alphainti*b[i][j][tt].rhor[k]+alphaintip*b[i][j][ll].rhol[k]);
                    // Calculation of total flux
                    if(vstar>0)
                    {
                        for(jj=0;jj<4;jj++)
                        {
                            b[i][j][ll].r[jj]=b[i][j][ll].r[jj]-Fcart[j][jj];
                            b[i][j][tt].r[jj]=b[i][j][tt].r[jj];
                        }
                        b[i][j][ll].r[1]=b[i][j][ll].r[1]-b[i][j][0].deltaeff[k]*pint*b[i][j][ll].nn[k].x*mod1(b[i][j][ll].area[k]);
                        b[i][j][ll].r[2]=b[i][j][ll].r[2]-b[i][j][0].deltaeff[k]*pint*b[i][j][ll].nn[k].y*mod1(b[i][j][ll].area[k]);
                    }
                    else
                    {
                        for(jj=0;jj<4;jj++)
                        {
                            b[i][j][ll].r[jj]=b[i][j][ll].r[jj];
                            b[i][j][tt].r[jj]=b[i][j][tt].r[jj]-Fcart[j][jj];
                        }
                        b[i][j][ll].r[1]=b[i][j][ll].r[1]-b[i][j][0].deltaeff[k]*pint*b[i][j][ll].nn[k].x*mod1(b[i][j][ll].area[k]);
                        b[i][j][ll].r[2]=b[i][j][ll].r[2]-b[i][j][0].deltaeff[k]*pint*b[i][j][ll].nn[k].y*mod1(b[i][j][ll].area[k]);
                    }
                }
            }
        }
    }
    for(int j=0;j<(ny-1);j++)
    {
        delete[] Flocal[j];
        delete[] Fcart[j];
    }
    delete Flocal;
    delete Fcart;
}
// Since we use the three stage TVD Runge Kutta method, within each time step, we have three intermediate steps which are represented by three solver functions
void solver1(cell***b,double deltat,int nx,int ny)
{
    int j;
    for(int l=0;l<2;l++)
    {
        #pragma omp parallel for private(j)
        for(int i=0;i<(nx-1);i++)
        {
            for(j=0;j<(ny-1);j++)
            {
                b[i][j][l].U[0]=b[i][j][l].U[0]+deltat*b[i][j][l].r[0]/b[i][j][l].volume;
                b[i][j][l].U[1]=b[i][j][l].U[1]+deltat*b[i][j][l].r[1]/b[i][j][l].volume;
                b[i][j][l].U[2]=b[i][j][l].U[2]+deltat*b[i][j][l].r[2]/b[i][j][l].volume;
                b[i][j][l].U[3]=b[i][j][l].U[3]+deltat*b[i][j][l].r[3]/b[i][j][l].volume;
            }
        }
    }
}
void solver2(cell***b,double deltat,int nx,int ny)
{
    int j;
    for(int l=0;l<2;l++)
    {
        #pragma omp parallel for private(j)
        for(int i=0;i<(nx-1);i++)
        {
            for(j=0;j<(ny-1);j++)
            {
                b[i][j][l].U[0]=0.75*b[i][j][l].W[0]+0.25*b[i][j][l].U[0]+0.25*(deltat/b[i][j][l].volume)*b[i][j][l].r[0];
                b[i][j][l].U[1]=0.75*b[i][j][l].W[1]+0.25*b[i][j][l].U[1]+0.25*(deltat/b[i][j][l].volume)*b[i][j][l].r[1];
                b[i][j][l].U[2]=0.75*b[i][j][l].W[2]+0.25*b[i][j][l].U[2]+0.25*(deltat/b[i][j][l].volume)*b[i][j][l].r[2];
                b[i][j][l].U[3]=0.75*b[i][j][l].W[3]+0.25*b[i][j][l].U[3]+0.25*(deltat/b[i][j][l].volume)*b[i][j][l].r[3];
            }
        }
    }
}
void solver3(cell***b,double deltat,int nx,int ny)
{
    // Here in this function,we finally get the values at the next time step. W stores the values at this time step which is then further used in intermiediate calcuations of the next time step
    int j;
    for(int l=0;l<2;l++)
    {
        #pragma omp parallel for private(j)
        for(int i=0;i<(nx-1);i++)
        {
            for (j=0;j<(ny-1);j++)
            {
                b[i][j][l].W[0]=b[i][j][l].U[0]=b[i][j][l].W[0]/3+2*b[i][j][l].U[0]/3+2*(deltat/b[i][j][l].volume)*b[i][j][l].r[0]/3;
                b[i][j][l].W[1]=b[i][j][l].U[1]=b[i][j][l].W[1]/3+2*b[i][j][l].U[1]/3+2*(deltat/b[i][j][l].volume)*b[i][j][l].r[1]/3;
                b[i][j][l].W[2]=b[i][j][l].U[2]=b[i][j][l].W[2]/3+2*b[i][j][l].U[2]/3+2*(deltat/b[i][j][l].volume)*b[i][j][l].r[2]/3;
                b[i][j][l].W[3]=b[i][j][l].U[3]=b[i][j][l].W[3]/3+2*b[i][j][l].U[3]/3+2*(deltat/b[i][j][l].volume)*b[i][j][l].r[3]/3;
            }
        }
    }
}

// This function calculates the values of primitive variables and other variables from the conservative variables at the end of each time step
void decode(cell***b,int nx,int ny,double emin,double radius,double eee)
{
    //double A[nx-1][ny-1][2],a[nx-1][ny-1][2],B[nx-1][ny-1],c[nx-1][ny-1];
    // Details of A, a, B, c may be found in Liou et al.
    double***A=new double**[nx-1];
    double***a=new double**[nx-1];
    double**B=new double*[nx-1];
    double**c=new double*[nx-1];
    double Adet,f1,f2,r,e2,gg;
    for(int i=0;i<(nx-1);i++)
    {
        A[i]=new double*[ny-1];
        a[i]=new double*[ny-1];
        B[i]=new double[ny-1];
        c[i]=new double[ny-1];
        for(int j=0;j<(ny-1);j++)
        {
            A[i][j]=new double[2];
            a[i][j]=new double[2];
        }
    }
    int j;
    //cout<<"a";
    // Details of the formulas may be found in the text
    for(int l=0;l<2;l++)
    {
        #pragma omp parallel for private(j)
        for(int i=0;i<(nx-1);i++)
        {
            for(j=0;j<(ny-1);j++)
            {
                A[i][j][l]=(b[i][j][l].U[3]-((pow(b[i][j][l].U[1],2)+pow(b[i][j][l].U[2],2))/(2*b[i][j][l].U[0])))*(b[i][j][l].gamma-1);
                a[i][j][l]=b[i][j][l].gamma*b[i][j][l].pinf+(b[i][j][l].gamma-1)*b[i][j][l].pint;
            }
        }
    }
    #pragma omp parallel for private(j)
    for(int i=0;i<(nx-1);i++)
    {
        for(j=0;j<(ny-1);j++)
        {
             B[i][j]=A[i][j][0]-a[i][j][0]+A[i][j][1]-a[i][j][1];
             c[i][j]=a[i][j][0]*A[i][j][1]+a[i][j][1]*A[i][j][0]-a[i][j][0]*a[i][j][1];
             //cout<<B[i][j]<<"\t"<<i<<"\t"<<j<<"\n";
        }
    }
    // Calculation of pressure and volume fraction of each phase takes place in this section
    #pragma omp parallel for private(j,Adet,f1,f2,r)
    for(int i=0;i<(nx-1);i++)
    {
        for(j=0;j<(ny-1);j++)
        {
            b[i][j][0].p=0.5*(B[i][j]+sqrt(pow(B[i][j],2)+4*c[i][j]));
            b[i][j][0].alpha=A[i][j][0]/(b[i][j][0].p+a[i][j][0]);
            //cout<<b[i][j][0].p<<"\n";
            // In order to improve accuracy of calculation of pressure and volume fraction, we employ Newton's iteration method for simultaneous
            // equations
            while(1)
                {
                    Adet=b[i][j][0].alpha*(a[i][j][0]-a[i][j][1])-b[i][j][0].p-a[i][j][0];  // Adet represents the determinant of the matrix used in Newton's Iteration method
                    f1=(b[i][j][0].p+a[i][j][0])*b[i][j][0].alpha-A[i][j][0];   // f1 and f2 are the values of first and second equation using the approximate values of p and volume fraction evaluated earlier
                    f2=(b[i][j][0].p+a[i][j][1])*(1-b[i][j][0].alpha)-A[i][j][1];
                    r=b[i][j][0].p;
                    // evaluation of p and alpha using newton's iteration method
                    b[i][j][0].p=b[i][j][0].p+((b[i][j][0].p+a[i][j][1])*f1+(b[i][j][0].p+a[i][j][0])*f2)/Adet;
                    b[i][j][0].alpha=b[i][j][0].alpha-((-1+b[i][j][0].alpha)*f1+b[i][j][0].alpha*f2)/Adet;
                    b[i][j][1].alpha=1-b[i][j][0].alpha;
                    b[i][j][1].p=b[i][j][0].p;
                    if(fabs(r-b[i][j][0].p)<pow(10,-5))
                        break;
                }
                if(b[i][j][0].alpha<emin)
                    b[i][j][0].alpha=emin;
                if(b[i][j][1].alpha<emin)
                    b[i][j][1].alpha=emin;

        }
    }
    // Calculation of remaining quantities by employing the equation of state and other equations
    for(int l=0;l<2;l++)
    {
        #pragma omp parallel for private(j)
        for(int i=0;i<(nx-1);i++)
        {
            for(j=0;j<(ny-1);j++)
            {
                b[i][j][l].rho=b[i][j][l].U[0]/b[i][j][l].alpha;
                b[i][j][l].vcart.x=b[i][j][l].U[1]/b[i][j][l].U[0];
                b[i][j][l].vcart.y=b[i][j][l].U[2]/b[i][j][l].U[0];
                b[i][j][l].E=(b[i][j][l].U[3]-b[i][j][l].pint*b[i][j][l].alpha)/b[i][j][l].U[0];
                b[i][j][l].T=(b[i][j][l].E-b[i][j][l].pinf/b[i][j][l].rho-pow(mod1(b[i][j][l].vcart),2)/2)*b[i][j][l].gamma/b[i][j][l].Cp;
                b[i][j][l].H=b[i][j][l].E+b[i][j][l].p/b[i][j][l].rho;
                b[i][j][l].a=sqrt(b[i][j][l].gamma*(b[i][j][l].p+b[i][j][l].pinf)/b[i][j][l].rho);
                b[i][j][l].e=(b[i][j][l].Cp*b[i][j][l].T/b[i][j][l].gamma)+b[i][j][l].pinf/b[i][j][l].rho;
                //cout<<b[i][j][l].rho<<"\n";
            }
        }
    }
    for(int i=0;i<(nx-1);i++)
    {
        for(int j=0;j<(ny-1);j++)
        {
            delete[] A[i][j];
            delete[] a[i][j];
        }
        delete[] A[i];
        delete[] a[i];
        delete[] B[i];
        delete[] c[i];
    }
    delete A;delete a;
    delete B;delete c;
}

// Whenever a particular phase is vanishing, we use this function. We implement this function whenever volume fraction of any phase
// becomes less than a specified tolerance value
void vanishing(cell***b,double emin,double emax,int nx,int ny)
{
    double e1,alpha1,g; int t;
    int j;
    for(int l=0;l<2;l++)
    {
        #pragma omp parallel for private(j,alpha1,e1,g,t)
        for(int i=0;i<(nx-1);i++)
        {
            for(j=0;j<(ny-1);j++)
            {
                // emax is the specified tolerance
                if(b[i][j][l].alpha<emax)
                {
                    if(l==1)
                        t=0;
                    else
                        t=1;
                    // if volume fraction becomes less than emin, then to preserve positivity we assign the value of emin to
                    // volume fraction
                    if(b[i][j][l].alpha<emin)
                    {
                        alpha1=emin;
                        b[i][j][l].alpha=emin;
                        e1=(alpha1-emin)/(emax-emin);
                        g=-pow(e1,2)*(2*e1-3);
                        // Here we blend the temperature and velocity of the vanishing phase with the dominant phase
                        b[i][j][l].vcart.x=g*b[i][j][l].vcart.x+(1-g)*b[i][j][t].vcart.x;
                        b[i][j][l].vcart.y=g*b[i][j][l].vcart.y+(1-g)*b[i][j][t].vcart.y;
                        b[i][j][l].T=g*b[i][j][l].T+(1-g)*b[i][j][t].T;
                    }
                    else
                    {
                        alpha1=b[i][j][l].alpha;
                        e1=(alpha1-emin)/(emax-emin);
                        g=-pow(e1,2)*(2*e1-3);
                        b[i][j][l].vcart.x=g*b[i][j][l].vcart.x+(1-g)*b[i][j][t].vcart.x;
                        b[i][j][l].vcart.y=g*b[i][j][l].vcart.y+(1-g)*b[i][j][t].vcart.y;
                        b[i][j][l].T=g*b[i][j][l].T+(1-g)*b[i][j][t].T;
                    }
                }
            }
        }
    }
}

// This function calculates the cfl number for each time step.
double stability(cell***b,double cfl,int nx,int ny)
{
    double deltat;
    double a=0.01;
    double tt;
    double c=0.01;
    int j;
    #pragma omp parallel for private(j,tt,c,a)
    for(int i=0;i<(nx-1);i++)
    {
        for(j=0;j<(ny-1);j++)
        {
            tt=b[i][j][0].widx/(max(b[i][j][0].a,b[i][j][1].a)+max(fabs(b[i][j][0].vcart.x),fabs(b[i][j][1].vcart.x)));
            a=min(a,tt);
            tt=b[i][j][0].widy/(max(b[i][j][0].a,b[i][j][1].a)+max(fabs(b[i][j][0].vcart.y),fabs(b[i][j][1].vcart.y)));
            c=min(c,tt);
        }
    }
    deltat=cfl*min(a,c);
    return deltat;
}

// This function prints output files after a certain number of time steps.
void writeoutput(cell***b,int nx,int ny,int gg,double ts)
{
    string filename = "ausmd2d1";
    string filetype = ".dat";
    stringstream ss;
    ofstream myfile;
    // Here we append the time at which we get an output file to the file name.
    ss<<filename<<ts<<filetype;
    myfile.open(ss.str().c_str());
    for(int j=1;j<(ny-2);j++)
    {
        for(int i=1;i<(nx-2);i++)
        {
            myfile<<b[i][j][0].center.x<<"\t"<<b[i][j][0].center.y<<"\t"<<b[i][j][0].rho<<"\t"<<b[i][j][0].p<<"\t"<<b[i][j][0].Tavg<<"\t"<<b[i][j][0].vavg.x<<"\t"<<b[i][j][0].vcart.y<<"\t"<<b[i][j][0].rhoavg<<"\t"<<b[i][j][0].alpha<<"\t"<<b[i][j][1].alpha<<"\t"<<b[i][j][1].rho<<"\t"<<mod1(b[i][j][0].delrho)<<"\n";
        }
    }
    myfile.close();
}
int main()
{
    int nx=901;    // specification of x nodes
    int ny=421;    // specification of y nodes
    double eee=pow(10,-3);
    vectorc**a=new vectorc*[nx];
    for(int i=0;i<nx;i++)
    {
        a[i]=new vectorc[ny];
    }
    // b represents the cell array which stores the values of variables and quantities for each cell
    cell***b=new cell**[nx-1];
    cell****bound1=new cell***[nx-1];
    cell****bound2=new cell***[nx-1];
    int gg=0;
    double Cpstar=2;          // Used in the calculation of interfacial pressure between different phases within each cell
    double emin=pow(10,-3);   // Specification of tolerance value for volume fraction
    double emax=pow(10,-1);
    double Kp=1;
    double Ku=1;
    //double pavgy[nx-1];
    //double uavgy[nx-1];
    //double rhoavgy[nx-1];
    //double perror[nx-1];
    //double rhoerror[nx-1];
    //double uerror[nx-1];
    double radius=3.2*pow(10,-3);     // radius of the bubble under consideration
    for(int i=0;i<(nx-1);i++)
    {
        b[i]=new cell*[ny-1];
        //cout<<i<<"\n";
        bound1[i]=new cell**[ny-1];
        bound2[i]=new cell**[ny-1];
        for(int j=0;j<(ny-1);j++)
        {
            bound1[i][j]=new cell*[2];
            bound2[i][j]=new cell*[2];
            for(int l=0;l<2;l++)
            {
                bound1[i][j][l]=new cell[4];
                bound2[i][j][l]=new cell[4];
            }
            b[i][j]=new cell[2];
        }
    }
    //cout<<"A"<<flush;
    int imin=0;
    int jmin=0;
    //double lengthx=10;
    //double lengthy=10;
     double deltax;
    double deltay;
    // Grid construction
    double llc=1.968*pow(10,-4);
    double llcp=1.8511*pow(10,-4);
    double cfl=0.2;
    int i;
    #pragma omp parallel for private(i,deltax)
    for(int j=0;j<ny;j++)
    {
        a[0][j].x=-15*pow(10,-3);
        for(i=1;i<nx;i++)
        {
            if(i<218)
            {
                deltax=(0.025-(i-215)*llc)*pow(10,-3);
                a[i][j].x=a[i-1][j].x+deltax;
            }
            else if(i>=218&&i<618)
            {
                deltax=0.025*pow(10,-3);
                a[i][j].x=a[i-1][j].x+deltax;
            }
            else
            {
                deltax=(0.025+(i-615)*llc)*pow(10,-3);
                a[i][j].x=a[i-1][j].x+deltax;
            }
        }
    }
    int j;
    #pragma omp parallel for private(j,deltay)
    for(i=0;i<nx;i++)
    {
        a[i][0].y=0;
        for(j=1;j<ny;j++)
        {
            if(j<200)
            {
                deltay=0.025*pow(10,-3);
                a[i][j].y=a[i][j-1].y+deltay;
            }
            else
            {
                deltay=(0.025+(j-198)*llcp)*pow(10,-3);
                a[i][j].y=a[i][j-1].y+deltay;
            }
        }
    }
    //cout<<"A"<<flush;
    for(int l=0;l<2;l++)
    {
        int t=0;int e=0;
        for(j=0;j<(ny-1);j++)
        {
            for(i=0;i<(nx-1);i++)
            {
                // specification of position vectors for vertices of each cell
                b[i][j][l].position[0]=a[t][e];
                b[i][j][l].position[1]=a[t+1][e];
                b[i][j][l].position[2]=a[t+1][e+1];
                b[i][j][l].position[3]=a[t][e+1];
                //cout<<b[i][j][l].position[0].x<<"\t"<<i<<"\t"<<j<<"\n";
                t=t+1;
                if(t==(nx-1))
                {
                    e=e+1;
                    t=0;
                }
            }
        }
    }
    for(i=0;i<nx;i++)
    {
        delete[] a[i];
    }
    delete a;
    geometry(b,nx,ny);
    //cout<<b[20][1][0].widx;
    // Specification of parameters for equation of state for each phase
    for(int l=0;l<2;l++)
    {
        #pragma omp parallel for private(j)
        for(i=0;i<(nx-1);j++)
        {
            for(j=0;j<(ny-1);j++)
            {
                if(l==0)
                {
                    b[i][j][l].gamma=1.4;
                    b[i][j][l].pinf=0;
                    b[i][j][l].Cp=1004.5;
                }
                else
                {
                    b[i][j][l].gamma=2.8;
                    b[i][j][l].pinf=8.5*pow(10,8);
                    b[i][j][l].Cp=4186;
                }
            }
        }
    }
    double e2, gg1;
    // Initialization of variables for each cell
    for(int l=0;l<2;l++)
    {
        #pragma omp parallel for private(j,e2,gg1)
        for(i=0;i<(nx-1);i++)
        {
            for(j=0;j<(ny-1);j++)
            {
                if(b[i][j][l].center.x<=-4*pow(10,-3))
                {
                    b[i][j][l].p=1.6*pow(10,9);
                    b[i][j][l].vcart.x=661.81;
                    b[i][j][l].vcart.y=0;
                    b[i][j][l].T=595.13;
                    if(l==0)
                        b[i][j][l].alpha=eee;
                    else
                        b[i][j][l].alpha=1-eee;
                }
                else if((pow(b[i][j][l].center.x,2)+pow(b[i][j][l].center.y,2))>pow(3.2*pow(10,-3),2))
                {
                    b[i][j][l].p=pow(10,5);
                    b[i][j][l].vcart.x=0;
                    b[i][j][l].vcart.y=0;
                    b[i][j][l].T=292.98;
                    if(l==0)
                        b[i][j][l].alpha=eee;
                    else
                        b[i][j][l].alpha=1-eee;
                }
                else
                {
                    b[i][j][l].p=pow(10,5);
                    b[i][j][l].vcart.x=0;
                    b[i][j][l].vcart.y=0;
                    b[i][j][l].T=292.98;
                    if(l==0)
                        b[i][j][l].alpha=1-eee;
                    else
                        b[i][j][l].alpha=eee;
                }
                b[i][j][l].rho=b[i][j][l].gamma*(b[i][j][l].p+b[i][j][l].pinf)/(b[i][j][l].Cp*b[i][j][l].T*(b[i][j][l].gamma-1));
                b[i][j][l].E=(b[i][j][l].Cp*b[i][j][l].T/b[i][j][l].gamma)+b[i][j][l].pinf/b[i][j][l].rho+pow(mod1(b[i][j][l].vcart),2)/2;
                b[i][j][l].H=b[i][j][l].E+b[i][j][l].p/b[i][j][l].rho;
                b[i][j][l].a=sqrt(b[i][j][l].gamma*(b[i][j][l].p+b[i][j][l].pinf)/b[i][j][l].rho);
                b[i][j][l].e=(b[i][j][l].Cp*b[i][j][l].T/b[i][j][l].gamma)+b[i][j][l].pinf/b[i][j][l].rho;
                if((sqrt(pow(b[i][j][0].center.x,2)+pow(b[i][j][0].center.y,2))>=(radius-2*0.025*pow(10,-3)))&&(sqrt(pow(b[i][j][0].center.x,2)+pow(b[i][j][0].center.y,2))<=(radius+2*0.025*pow(10,-3))))
                {
                    e2=(sqrt(pow(b[i][j][0].center.x,2)+pow(b[i][j][0].center.y,2))-(radius-2*0.025*pow(10,-3)))/(4*0.025*pow(10,-3));
                    gg1=-pow(e2,2)*(2*e2-3);
                    b[i][j][0].alpha=(gg1*eee+(1-gg1)*(1-eee));
                    b[i][j][1].alpha=1-b[i][j][0].alpha;
                }
                //cout<<b[0][0][0].a<<"\n";
            }
        }
    }
    double ts=0;int bb,cc,dd;bb=1;
    double deltat;
    // the main loop starts here
    for(;;)
    {
        if(bb==1)
        {
            intpressure(b,Cpstar,nx,ny);
            variable_calculator(b,nx,ny);
        }
        boundary(b,bound1,bound2,nx,ny,eee);
        localtocart(b,bound1,bound2,nx,ny);
        reconstruction(b,bound1,bound2,gg,nx,ny);
        intsoundspeed(b,nx,ny,eee);
        //cout<<b[250][0][0].rho<<"\n";
        interfacial(b,nx,ny,Kp,Ku,eee);
        //cout<<bound2[0][0][3].widx<<"\n";
        //if(ts<6.9*pow(10,-6))
        {
            deltat=3.125*pow(10,-10);
        }
        //else
        {
        //    deltat=stability(b,cfl,nx,ny);
        }
        if(bb==1)
        {
            solver1(b,deltat,nx,ny);
        }
        else if(cc==1)
        {
            solver2(b,deltat,nx,ny);
        }
        else if(dd==1)
        {
            solver3(b,deltat,nx,ny);
            gg=gg+1;
            ts=ts+deltat;
        }
        //cout<<"a";
        decode(b,nx,ny,emin,radius,eee);
        cout<<"A"<<flush;
        vanishing(b,emin,emax,nx,ny);

        if(bb==1)
        {
            bb=0;cc=1;dd=0;
        }
        else if(cc==1)
        {
            bb=0;cc=0;dd=1;
        }
        else if(dd==1)
        {
            bb=1;cc=0;dd=0;
        }
        //if(gg%100==0)
        if(((gg%400)==0)&&(bb==1))
        {
            for(int l=0;l<2;l++)
            {
                #pragma omp parallel for private(j)
                for(i=0;i<(nx-1);i++)
                {
                    for(j=0;j<(ny-1);j++)
                    {
                       b[i][j][l].Tavg=(b[i][j][0].alpha*b[i][j][0].rho*b[i][j][0].T+b[i][j][1].alpha*b[i][j][1].rho*b[i][j][1].T)/(b[i][j][0].alpha*b[i][j][0].rho+b[i][j][1].alpha*b[i][j][1].rho);
                       b[i][j][l].vavg.x=(b[i][j][0].alpha*b[i][j][0].rho*b[i][j][0].vcart.x+b[i][j][1].alpha*b[i][j][1].rho*b[i][j][1].vcart.x)/(b[i][j][0].alpha*b[i][j][0].rho+b[i][j][1].alpha*b[i][j][1].rho);
                       b[i][j][l].rhoavg=(b[i][j][0].alpha*b[i][j][0].rho+b[i][j][1].alpha*b[i][j][1].rho)/(b[i][j][0].alpha+b[i][j][1].alpha);
                    }
                }
            }
            #pragma omp parallel for private(j)
            for(i=1;i<(nx-2);i++)
            {
                for(j=1;j<(ny-2);j++)
                {
                    b[i][j][0].delrho.x=b[i][j][1].delrho.x=(((b[i+1][j][0].rhoavg*b[i+1][j][0].widx+b[i][j][0].rhoavg*b[i][j][0].widx)/(b[i+1][j][0].widx+b[i][j][0].widx))-((b[i][j][0].rhoavg*b[i][j][0].widx+b[i-1][j][0].rhoavg*b[i-1][j][0].widx)/(b[i][j][0].widx+b[i-1][j][0].widx)))/(b[i][j][0].widx);
                    b[i][j][0].delrho.y=b[i][j][1].delrho.y=(((b[i][j+1][0].rhoavg*b[i][j+1][0].widy+b[i][j][0].rhoavg*b[i][j][0].widy)/(b[i][j+1][0].widy+b[i][j][0].widy))-((b[i][j][0].rhoavg*b[i][j][0].widy+b[i][j-1][0].rhoavg*b[i][j-1][0].widy)/(b[i][j][0].widy+b[i][j-1][0].widy)))/(b[i][j][0].widy);
                }
            }
            writeoutput(b,nx,ny,gg,ts);
        }
          cout<<(ts/(5*pow(10,-6)))*100<<"\n";
        if(ts>5*pow(10,-6))
            break;
    }
    for(int l=0;l<2;l++)
    {
        #pragma omp parallel for private(j)
        for(i=0;i<(nx-1);i++)
        {
            for(j=0;j<(ny-1);j++)
            {
               b[i][j][l].Tavg=(b[i][j][0].alpha*b[i][j][0].rho*b[i][j][0].T+b[i][j][1].alpha*b[i][j][1].rho*b[i][j][1].T)/(b[i][j][0].alpha*b[i][j][0].rho+b[i][j][1].alpha*b[i][j][1].rho);
               b[i][j][l].vavg.x=(b[i][j][0].alpha*b[i][j][0].rho*b[i][j][0].vcart.x+b[i][j][1].alpha*b[i][j][1].rho*b[i][j][1].vcart.x)/(b[i][j][0].alpha*b[i][j][0].rho+b[i][j][1].alpha*b[i][j][1].rho);
               b[i][j][l].rhoavg=(b[i][j][0].alpha*b[i][j][0].rho+b[i][j][1].alpha*b[i][j][1].rho)/(b[i][j][0].alpha+b[i][j][1].alpha);
            }
        }
    }
    #pragma omp parallel for private(j)
    for(i=1;i<(nx-2);i++)
    {
        for(j=1;j<(ny-2);j++)
        {
            b[i][j][0].delrho.x=b[i][j][1].delrho.x=((b[i+1][j][0].rhoavg*b[i+1][j][0].widx+b[i][j][0].rhoavg*b[i][j][0].widx)/(b[i+1][j][0].widx+b[i][j][0].widx))-((b[i][j][0].rhoavg*b[i][j][0].widx+b[i-1][j][0].rhoavg*b[i-1][j][0].widx)/(b[i][j][0].widx+b[i-1][j][0].widx));
            b[i][j][0].delrho.y=b[i][j][1].delrho.y=((b[i][j+1][0].rhoavg*b[i][j+1][0].widy+b[i][j][0].rhoavg*b[i][j][0].widy)/(b[i][j+1][0].widy+b[i][j][0].widy))-((b[i][j][0].rhoavg*b[i][j][0].widy+b[i][j-1][0].rhoavg*b[i][j-1][0].widy)/(b[i][j][0].widy+b[i][j-1][0].widy));
        }
    }
    /*for(int i=0;i<(nx-1);i++)
    {
        pavgy[i]=0;rhoavgy[i]=0;uavgy[i]=0;
        for(int j=0;j<(ny-1);j++)
        {
            pavgy[i]=pavgy[i]+b[i][j][0].p*b[i][j][0].widy/lengthy;
            rhoavgy[i]=rhoavgy[i]+b[i][j][0].rhoavg*b[i][j][0].widy/lengthy;
            uavgy[i]=uavgy[i]+b[i][j][0].vavg.x*b[i][j][0].widy/lengthy;
        }
    }
    double sump;double sumrho;double sumu;
    for(int i=0;i<(nx-1);i++)
    {
        sump=0;sumrho=0;sumu=0;
        for(int j=0;j<(ny-1);j++)
        {
            sump=sump+(pow(b[i][j][0].p-pavgy[i],2))*b[i][j][0].widy;
            sumrho=sumrho+(pow(b[i][j][0].rhoavg-rhoavgy[i],2))*b[i][j][0].widy;
            sumu=sumu+(pow(b[i][j][0].vavg.x-uavgy[i],2))*b[i][j][0].widy;
        }
        perror[i]=sqrt(sump/lengthy);
        rhoerror[i]=sqrt(sumrho/lengthy);
        uerror[i]=sqrt(sumu/lengthy);
    }*/
    writeoutput(b,nx,ny,gg,ts);
    for(i=0;i<(nx-1);i++)
    {
        for(j=0;j<(ny-1);j++)
        {
            for(int l=0;l<2;l++)
            {
                delete[] bound1[i][j][l];
                delete[] bound2[i][j][l];
            }
            delete[] b[i][j];
            delete[] bound1[i][j];
            delete[] bound2[i][j];
        }
        delete[] bound1[i];
        delete[] bound2[i];
        delete[] b[i];
    }
    delete b;
    delete bound1;
    delete bound2;
    return 0;
}
