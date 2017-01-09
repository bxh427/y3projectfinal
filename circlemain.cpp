/* Ben Handley - 1423327
Mathematical Billiards
Circular Geometry
This code produces trajectory data for a billiard on a circular table.
The coordinate axes has origin at centre of circle. */

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

int main()
{
	// Both angles measured in multiples of pi
	double theta;			// Starting bounce position measured as an angle from the positive x axis
	double alpha;			// Bounce angle measured from the tangent of the circle
	double bounces;			// Total number of bounces to be taken
	double pi = acos(-1);	// Define pi
	
	// Define file streams for user and gnuplot
	ofstream outs;
	ofstream gnu;
	outs.open("Circle.txt");
	gnu.open("gnuCircle.txt");
	
	// User interface to obtain starting conditions and number of bounces
	cout << "This program calculates subsequent positions of bounces of a billiard ball within a circular table."
		 << "\n" << "The circular table is set at the origin and the results are independent of the table radius."
		 << "\n" << "The angle theta is the angle measured from the positive x axis to the ball's location in radians."
		 << "\n" << "The angle alpha is the angle measured from the tangent of the circle and indicates the"
		 << "\n" << "direction of the ball's bounce. It is measured anticlockwise from the tangent in radians."
		 << "\n" << "Please input the initial value of theta in units of pi and press [rtn]: ";
	cin  >> theta;
	cout << "Please input the value of alpha in units of pi and press [rtn]: ";
	cin  >> alpha;
	cout << "Please input the number of bounces required and press [rtn]: ";
	cin	 >> bounces;
	
	// Output headers
	cout << "Calculating, please wait..." << endl;
	
	outs << "\n" << "For a billiard ball starting at a position "
		 << theta << " pi from the positive x axis and making bounces of angle "
		 << alpha << " pi, the following path is taken: \n"
		 << setw(10) << "Bounce"
		 << setw(25) << "x position"
		 << setw(25) << "y position" << endl;
	
	gnu  << setw(25) << "#x"
		 << setw(25) << "y" << endl;

	// For loop to iterate over each bounce
	for(int n=0; n<=bounces; n++)
	{
		// Output results	 
		outs << setprecision(18) << setw(10) << n
			 << setw(25) << cos(theta*pi)
			 << setw(25) << sin(theta*pi) << endl;
		
		gnu  << setprecision(18) << setw(25) << cos(theta*pi)
			 << setw(25) << sin(theta*pi) << endl;
		
		// Move theta on by two alpha to produce next position
		theta+=2*alpha;
		
		// Enforce periodicity
		if(theta>=2) theta-=2;
	}
	
	cout << "\nResults complete.\n"
		 << "To view the results, open 'Circle.txt'. \n"
		 << "To plot this data, load 'circplot.gnu'.\n";
	
	outs.close();
	gnu.close();
	
	return 0;
}

