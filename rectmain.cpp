/* Ben Handley - 1423327
Mathematical Billiards
Rectangular Geometry
This code produces trajectory data for a billiard on a rectangular table of user-defined size.
The coordinate axes has origin at bottom left corner of rectangle. */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

// Define a vector struct
struct Vector
{
	double xcomp;
	double ycomp;
};

// Prototyping functions
Vector getInitPos(char SIDEINIT,double A,double B);
Vector getNorm(Vector POS,double A,double B);
double DotProduct(Vector vec1,Vector vec2);
Vector ScalVec(double scal,Vector vec);
Vector VecOut(Vector VecIn,Vector Vec2);
Vector TIncSide(Vector Vinc,Vector Pos,double A,double B);
Vector BIncSide(Vector Vinc,Vector Pos,double A,double B);
Vector LIncSide(Vector Vinc,Vector Pos,double A,double B);
Vector RIncSide(Vector Vinc,Vector Pos,double A,double B);

int main()
{
	double a,b;		// Height and width of rectangle
	char sideinit;	// Side from which ball starts
	Vector InitPos;	// Initial position
	Vector V_in;	// Incoming vector
	Vector Norm;	// Normal at surface
	double VdotN;	// Vector dotted with normal
	Vector VNN;		// Dot product multiplied by normal
	Vector V_out;	// Outgoing vector
	Vector IncPos;	// Incident position
	int bounces;	// Total number of bounces to be taken
	
	// Define file streams for user and gnuplot
	ofstream outs;
	ofstream gnu;
	outs.open("Rectangle.txt");
	gnu.open("gnuRectangle.txt");
	
	// User interface to obtain starting conditions and number of bounces
	cout << "Input the height of the rectangle: ";
	cin  >> a;
	cout << "Input the width of the rectangle: ";
	cin  >> b;
	cout << "Input which side the ball starts from according to:"
		 << "\n Top = t"
		 << "\n Bottom = b"
		 << "\n Left = l"
		 << "\n Right = r" << endl;
	cin  >> sideinit;
	InitPos = getInitPos(sideinit,a,b);
	cout << "Input the x component of the incoming vector: ";
	cin  >> V_in.xcomp;
	cout << "Input the y component of the incoming vector: ";
	cin  >> V_in.ycomp;
	cout << "Input the number of bounces required: ";
	cin  >> bounces;

	// Output headers and initial conditions
	cout << "Calculating, please wait..." << endl;
		 
	outs << "\n" << "For a billiard ball starting at position ("
		 << InitPos.xcomp << "," << InitPos.ycomp << ") and with initial vector ("
		 << V_in.xcomp << "," << V_in.ycomp << ") the following path is taken: \n"
		 << setw(6)  << "Bounce"
		 << setw(20) << "x"
		 << setw(20) << "y" << endl
		 << setw(6) << 0
		 << setw(20) << InitPos.xcomp
		 << setw(20) << InitPos.ycomp << endl;
		 
	gnu  << setw(25) << "#x"
		 << setw(25) << "y" << endl
		 << setw(25) << InitPos.xcomp
		 << setw(25) << InitPos.ycomp << endl;
	
	// For loop to iterate over each bounce
	for(int n=1;n<=bounces;n++)
	{
		// Find normal at surface
		Norm = getNorm(InitPos,a,b);
		
		// Dot incoming vector with normal
		VdotN = DotProduct(V_in,Norm);
		
		// Multiply dot product by normal
		VNN = ScalVec(VdotN,Norm);
		
		// Find outward vector
		V_out = VecOut(V_in,VNN);
		
		// Find the incident point by running function appropriate to current side
		if(InitPos.ycomp==a)
		{
			IncPos = TIncSide(V_out,InitPos,a,b);
		}
		
		else if(InitPos.ycomp==0)
		{
			IncPos = BIncSide(V_out,InitPos,a,b);
		}
		
		else if(InitPos.xcomp==0)
		{
			IncPos = LIncSide(V_out,InitPos,a,b);
		}
		
		else if(InitPos.xcomp==b)
		{
			IncPos = RIncSide(V_out,InitPos,a,b);
		}
		
		// Print results
		outs << setw(6) << n
			 << setw(20) << IncPos.xcomp
			 << setw(20) << IncPos.ycomp << endl;
		gnu  << setw(25) << IncPos.xcomp
			 << setw(25) << IncPos.ycomp << endl;
			 
		// Reset conditions	 
		InitPos = IncPos;
		V_in = V_out;
	}
	
	cout << "\nResults complete.\n"
		 << "To view this data, open 'Rectangle.txt'.\n"
		 << "To plot this data, open 'rectplot.gnu' and set a and b equal to the same height and width respectively."
		 << "\nThen load this file in gnuplot.\n";
	
	outs.close();
	gnu.close();
	
	return 0;
}

// Function to determine initial position depending on specified side
Vector getInitPos(char SIDEINIT, double A, double B)
{
	Vector INITPOS;

	// If initial side is the top
	if(SIDEINIT=='t')
	{
		// Set y component to be the maximum
		INITPOS.ycomp = A;
		
		// Request x coordinate
		cout << "Input the x coordinate between 0 and " << B << " : ";
		cin  >> INITPOS.xcomp;
	}
	
	// If initial side is the bottom
	if(SIDEINIT=='b')
	{
		// Set y component to be zero
		INITPOS.ycomp = 0;
		
		// Request x coordinate
		cout << "Input the x coordinate between 0 and " << B << " : ";
		cin  >> INITPOS.xcomp;
	}

	// If initial side is the left
	if(SIDEINIT=='l')
	{
		// Set x component to be zero
		INITPOS.xcomp = 0;
		
		// Request y coordinate
		cout << "Input the y coordinate between 0 and " << A << " : ";
		cin  >> INITPOS.ycomp;
	}

	// If initial side is the right
	if(SIDEINIT=='r')
	{
		// Set x component to be the maximum
		INITPOS.xcomp = B;
		
		// Request y coordinate
		cout << "Input the y coordinate between 0 and " << A << " : ";
		cin  >> INITPOS.xcomp;
	}
	
	return INITPOS;
}

// Function to return normal depending on current side (there are 4 discrete choices of normal)
Vector getNorm(Vector POS, double A, double B)
{
	Vector NORM;
	if(POS.xcomp==0)
	{
		NORM.xcomp = 1;
		NORM.ycomp = 0;
	}
	
	else if(POS.xcomp==B)
	{
		NORM.xcomp = -1;
		NORM.ycomp = 0;
	}
	
	else if(POS.ycomp==0)
	{
		NORM.xcomp = 0;
		NORM.ycomp = 1;
	}
	
	else if(POS.ycomp==A)
	{
		NORM.xcomp = 0;
		NORM.ycomp = -1;
	}
	
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

// Function to determine intercept if starting on top side
Vector TIncSide(Vector Vinc,Vector Pos,double A,double B)
{
	Vector INTERCEPT;
	double grad = Vinc.ycomp/Vinc.xcomp;
	double L;
	
	// Left hand intercept
	if(grad > 0 && grad < (A/Pos.xcomp))
	{
		INTERCEPT.xcomp = 0;
		L = -1*(Pos.xcomp/Vinc.xcomp);
		INTERCEPT.ycomp = Pos.ycomp + L*Vinc.ycomp;
	}
	
	// Right hand intercept
	else if(grad < 0 && grad > (A/(Pos.xcomp-B)))
	{
		INTERCEPT.xcomp = B;
		L = (B-Pos.xcomp)/Vinc.xcomp;
		INTERCEPT.ycomp = Pos.ycomp + L*Vinc.ycomp;
	}
	
	// Bottom intercept
	else
	{
		INTERCEPT.ycomp = 0;
		L = -1*(Pos.ycomp/Vinc.ycomp);
		INTERCEPT.xcomp = Pos.xcomp + L*Vinc.xcomp;
	}
	
	return INTERCEPT;
}

// Function to determine intercept if starting on bottom side
Vector BIncSide(Vector Vinc,Vector Pos,double A,double B)
{
	Vector INTERCEPT;
	double grad = Vinc.ycomp/Vinc.xcomp;
	double L;
	
	// Right hand intercept
	if(grad > 0 && grad < (A/(B-Pos.xcomp)))
	{
		INTERCEPT.xcomp = B;
		L = (B-Pos.xcomp)/Vinc.xcomp;
		INTERCEPT.ycomp = Pos.ycomp + L*Vinc.ycomp;
	}
	
	// Left hand intercept
	else if(grad < 0 && grad > (-1*A/Pos.xcomp))
	{
		INTERCEPT.xcomp = 0;
		L = -1*(Pos.xcomp/Vinc.xcomp);
		INTERCEPT.ycomp = Pos.ycomp + L*Vinc.ycomp;
	}
	
	// Top intercept
	else
	{
		INTERCEPT.ycomp = A;
		L = (A-Pos.ycomp)/Vinc.ycomp;
		INTERCEPT.xcomp = Pos.xcomp + L*Vinc.xcomp;
	}
	
	return INTERCEPT;
}

// Function to determine intercept if starting on left side
Vector LIncSide(Vector Vinc,Vector Pos,double A,double B)
{
	Vector INTERCEPT;
	double grad = Vinc.ycomp/Vinc.xcomp;
	double L;
	
	// Top intercept
	if(grad > ((A-Pos.ycomp)/B))
	{
		INTERCEPT.ycomp = A;
		L = (A-Pos.ycomp)/Vinc.ycomp;
		INTERCEPT.xcomp = Pos.xcomp + L*Vinc.xcomp;
	}
	
	// Bottom intercept
	else if(grad < (-1*Pos.ycomp/B))
	{
		INTERCEPT.ycomp = 0;
		L = -1*(Pos.ycomp/Vinc.ycomp);
		INTERCEPT.xcomp = Pos.xcomp + L*Vinc.xcomp;
	}
	
	// Right hand intercept
	else
	{
		INTERCEPT.xcomp = B;
		L = (B-Pos.xcomp)/Vinc.xcomp;
		INTERCEPT.ycomp = Pos.ycomp + L*Vinc.ycomp;
	}
	
	return INTERCEPT;
}

// Function to determine intercept if starting on right side
Vector RIncSide(Vector Vinc,Vector Pos,double A,double B)
{
	Vector INTERCEPT;
	double grad = Vinc.ycomp/Vinc.xcomp;
	double L;
	
	// Bottom intercept
	if(grad > (Pos.ycomp/B))
	{
		INTERCEPT.ycomp = 0;
		L = -1*(Pos.ycomp/Vinc.ycomp);
		INTERCEPT.xcomp = Pos.xcomp + L*Vinc.xcomp;
	}
	
	// Top intercept
	else if(grad < 0 && grad < ((Pos.ycomp-A)/B))
	{
		INTERCEPT.ycomp = A;
		L = (A-Pos.ycomp)/Vinc.ycomp;
		INTERCEPT.xcomp = Pos.xcomp + L*Vinc.xcomp;
	}
	
	// Left intercept
	else
	{
		INTERCEPT.xcomp = 0;
		L = -1*(Pos.xcomp/Vinc.xcomp);
		INTERCEPT.ycomp = Pos.ycomp + L*Vinc.ycomp;
	}
	
	return INTERCEPT;
}
