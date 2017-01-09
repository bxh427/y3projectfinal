/* Ben Handley - 1423327
Mathematical Billiards
Circular Geometry
This code produces Poincare data for a billiard on a circular table.
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
	
	// Define file streams for user and gnuplot
	ofstream gnu;
	ofstream outs;
	outs.open("poincareCircle.txt");
	gnu.open("gnupoincareCircle.txt");
	
	// User interface to obtain starting conditions and number of bounces
	cout << "This program calculates subsequent positions of bounces of a billiard ball within a circular table."
		 << "\n" << "The circular table is set at the origin and the results are independent of the table radius."
		 << "\n" << "The angle theta is the angle measured from the positive x axis to the ball's location in radians."
		 << "\n" << "The angle alpha is the angle measured from the tangent of the circle and indicates the"
		 << "\n" << "direction of the ball's bounce. It is measured anticlockwise from the tangent in radians."
		 << "\n" << "Please input the initial value of theta in units of pi and press [rtn]: ";
	cin  >> theta;
	cout << "Please input the number of bounces required and press [rtn]: ";
	cin	 >> bounces;
	
	// Initial value of alpha
	alpha = 0;

	// Output headers
	cout << "Calculating, please wait..." << endl;
	outs << "\n" << "For a billiard ball starting at a position "
		 << theta << " pi from the positive x axis the following data is produced for various alpha: " << endl
		 << setw(10) << "Bounce"
		 << setw(25) << "theta"
		 << setw(25) << "alpha" << endl;
	
	gnu  << setw(25) << "#theta"
		 << setw(25) << "alpha" << endl;
	
	// Define a "permanent" theta which the system can be reset to
	double permtheta = theta;

	// For loop to iterate over each alpha
	for(int i=0; i<=10; i++)
	{
		// For loop to iterate over each bounce
		for(int n=0; n<=bounces; n++)
		{
			// Output results
			outs << setprecision(15) << setw(10) << n
				 << setw(25) << theta
				 << setw(25) << alpha << endl;
		
			gnu  << setprecision(15) << setw(25) << theta
				 << setw(25) << alpha << endl;
		
			// Move theta on by two alpha to produce next position
			if(alpha >= 0) theta+=2*alpha;
			if(alpha < 0) theta-=2*alpha; 		

			// Enforce periodicity
			while(theta>=2)
			{
				theta-=2;
			}
		}
		
		// Reset initial position
		theta = permtheta;

		// If conditions to set the next alpha
		if(alpha >= 0 && alpha < 0.5)
		{
			alpha += 0.1;
		}

		else if(alpha >= 0.5)
		{
			alpha = -0.1;
		}

		else if(alpha < 0 && alpha >= -0.5)
		{
			alpha -= 0.1;
		}

		else if(alpha < -0.5)
		{
			alpha = 0.1;
		}
	}
	
	cout << "\nResults complete.\n"
		 << "To view this data, open 'poincareCircle.txt'. \n"
		 << "To plot this data, load 'circpoinplot.gnu'.\n";
	
	gnu.close();
	outs.close();
	
	return 0;
}

