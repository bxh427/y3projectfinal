/* Ben Handley - 1423327
Mathematical Billiards
Elliptical Geometry
This code produces trajectory data for a billiard on a elliptical table of user-defined size.
The coordinate axes has origin at centre of ellipse. */

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

// Define a vector struct
struct Vector
{
	double xcomp;
	double ycomp;
};

// Prototyping functions
double Radius(double major, double minor, double angle);
Vector Normal(double X_reduced, double Y_reduced);
double DotProduct(Vector vec1,Vector vec2);
Vector ScalVec(double scal,Vector vec);
Vector VecOut(Vector VecIn,Vector Vec2);
double calcLamb(double major, double minor, Vector vec, Vector pos);
Vector IncPos(Vector initpos, double L, Vector vout);

int main()
{
	double thetainit;		// Starting bounce position measured as an angle from the positive x axis
	int bounces;			// Total number of bounces to be taken
	double pi = acos(-1);	// Define pi
	double a,b;				// Semimajor and semiminor axes
	Vector V_in;			// Incoming vector
	Vector Pos1, redPos1;	// Initial position and "reduced" position	
	Vector Norm;			// Normal at surface
	double VdotN;			// Vector dotted with normal
	Vector VNN;				// Dot product multiplied by normal
	Vector V_out;			// Outgoing vector
	double Lambda;			// Lambda in vector equation of a line
	Vector Pos2;			// New position
	
	// Define file streams for user and gnuplot
	ofstream outs;
	ofstream gnu;
	outs.open("Ellipse.txt");
	gnu.open("gnuEllipse.txt");
	
	// User interface to obtain starting conditions and number of bounces
	cout << "This program calculates subsequent positions of bounces of a billiard ball within an elliptical table."
		 << "\n" << "The elliptical table is set at the origin."
		 << "\n" << "The angle theta is the angle measured from the positive x axis to the ball's location in radians."
		 << "\n" << "Please input the initial value of theta in units of pi and press [rtn]: ";
	cin >> thetainit;
	cout << "Input the semimajor axis a: ";
	cin >> a;
	cout << "Input the semiminor axis b: ";
	cin >> b;
	cout << "Input the x component of the initial vector direction: ";
	cin >> V_in.xcomp;
	cout << "Input the y component of the initial vector direction: ";
	cin >> V_in.ycomp;
	cout << "Input the number of bounces required: ";
	cin >> bounces;

	// Finding initial x and y for given initial conditions
	Pos1.xcomp = Radius(a,b,thetainit*pi) * cos(thetainit*pi);
	Pos1.ycomp = Radius(a,b,thetainit*pi) * sin(thetainit*pi);
	
	// Output headers and initial conditions
	cout << "Calculating, please wait..." << endl;
		 
	outs << "\n" << "For a billiard ball starting at a position "
		 << thetainit << " pi from the positive x axis and with initial vector ("
		 << V_in.xcomp << "," << V_in.ycomp << ") the following path is taken: \n"
		 << setw(10) << "Bounce"
		 << setw(25) << "x"
		 << setw(25) << "y" << endl
		 << setw(10) << 0
		 << setw(25) << Pos1.xcomp
		 << setw(25) << Pos1.ycomp << endl;
	
	gnu	 << setw(25) << "#x"
		 << setw(25) << "y" << endl
		 << setw(25) << Pos1.xcomp
		 << setw(25) << Pos1.ycomp << endl;
		
	// For loop to iterate over each bounce	
	for(int n=1; n<=bounces; n++)
	{
		// Define "reduced" variables Xinit, Yinit
		redPos1.xcomp = Pos1.xcomp/(a*a);
		redPos1.ycomp = Pos1.ycomp/(b*b);
		
		// Define vector normal to ellipse surface
		Norm = Normal(redPos1.xcomp, redPos1.ycomp);
				
		// Dot product of vector with normal
		VdotN = DotProduct(V_in, Norm);
		
		// Multiply dot product with normal
		VNN = ScalVec(VdotN,Norm);
		
		// Find reflected vector
		V_out = VecOut(V_in,VNN);
		
		// Calculate lambda for vector equation of a line
		Lambda = calcLamb(a,b,V_out,Pos1);
		
		// Find bounce point
		Pos2 = IncPos(Pos1, Lambda, V_out);
		
		// Print results
		outs << setprecision(18) << setw(10) << n
			 << setw(25) << Pos2.xcomp
			 << setw(25) << Pos2.ycomp << endl;
		
		gnu  << setprecision(18)
			 << setw(25) << Pos2.xcomp
			 << setw(25) << Pos2.ycomp << endl;
		
		// Reset conditions
		Pos1 = Pos2;
		V_in = V_out;
	}
	
	cout << "\nResults complete.\n"
		 << "To view this data, open 'Ellipse.txt'.\n"
		 << "To plot this data, open 'ellplot.gnu' and set a and b equal to the same semi-major and -minor axes respectively."
		 << "\nThen load this file in gnuplot.\n";
	
	outs.close();
	gnu.close();
	
	return 0;
}

// Function which returns the radius from the origin for a given angle
double Radius(double major, double minor, double angle)
{
	return (minor*major)/sqrt((minor*minor*cos(angle)*cos(angle))+(major*major*sin(angle)*sin(angle)));
}

// Function for calculating normal
Vector Normal(double X_reduced, double Y_reduced)
{
	Vector NORM;
	NORM.xcomp = -X_reduced/sqrt(X_reduced*X_reduced + Y_reduced*Y_reduced);
	NORM.ycomp = -Y_reduced/sqrt(X_reduced*X_reduced + Y_reduced*Y_reduced);
	
	return NORM;
}

// General function for producing a dot product of two vectors
double DotProduct(Vector vec1, Vector vec2)
{
	return (vec1.xcomp*vec2.xcomp) + (vec1.ycomp*vec2.ycomp);
}

// General function for multiplying a vector by a scalar
Vector ScalVec(double scal, Vector vec)
{
	Vector scalvec;
	scalvec.xcomp = scal*vec.xcomp;
	scalvec.ycomp = scal*vec.ycomp;
	
	return scalvec;
}

// Function to calculate outward vector
Vector VecOut(Vector VecIn, Vector Vec2)
{
	Vector VOUT;
	VOUT.xcomp = VecIn.xcomp - 2*Vec2.xcomp;
	VOUT.ycomp = VecIn.ycomp - 2*Vec2.ycomp;
	
	return VOUT;
}

// Function to calculate lambda for vector equation of a line
double calcLamb(double major, double minor, Vector vec, Vector pos)
{
	return -2*((minor*minor*vec.xcomp*pos.xcomp + major*major*vec.ycomp*pos.ycomp)/(vec.xcomp*vec.xcomp*minor*minor + vec.ycomp*vec.ycomp*major*major));
}

// Function to find new incident position
Vector IncPos(Vector initpos, double L, Vector vout)
{
	Vector POS;
	POS.xcomp = initpos.xcomp + L*vout.xcomp;
	POS.ycomp = initpos.ycomp + L*vout.ycomp;
		
	return POS;
}
