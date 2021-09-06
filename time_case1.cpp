#include<bits/stdc++.h>
#define e 2.718281828459
using namespace std;

double df(double y, int a)//This function is useful in shortening of Adomian Polynomials
{
	double b=1;
	for(int i=0;i<a;i++)
	{
		b*=(y-i);
	}
	return b;
}

double masterpoly(double x, double u[], int order)//This function returns the value of Adomian polynomial
{
	double polymaster[10];
	if(order == 0)
	{
		polymaster[0] = pow(u[0],x);
		return polymaster[0];
	}
	else if(order == 1)
	{	
		polymaster[1] = x*pow(u[0],x-1)*u[1];
		return polymaster[1];
	}
	else if(order == 2)
	{
		polymaster[2] = x*pow(u[0],x-1)*u[2] + (0.5*x*(x-1)*pow(u[0],x-2)*u[1]*u[1]);
		return polymaster[2];
	}
	else if(order == 3)
	{
		polymaster[3] = x*pow(u[0],x-1)*u[3] + (x*(x-1)*pow(u[0],x-2)*u[1]*u[2]) + ((df(x,3)*pow(u[0],x-3)*pow(u[1],3))/6);
		return polymaster[3];		
	}
	else if(order == 4)
	{
		polymaster[4] = x*pow(u[0],x-1)*u[4] + (x*(x-1)*pow(u[0],x-2)*((u[2]*u[2]/2)+(u[1]*u[3]))) + ((x*(x-1)*(x-2)*pow(u[0],x-3)*pow(u[1],2)*u[2])/2) + ((x*(x-1)*(x-2)*(x-3)*pow(u[0],x-4)*pow(u[1],4))/24);
		return polymaster[4];
	}
	else if(order == 5)	
	{
		polymaster[5] = x*pow(u[0],x-1)*u[5] + (x*(x-1)*pow(u[0],x-2)*((u[2]*u[3])+(u[1]*u[4]))) + ((df(x,3)*pow(u[0],x-3)*( (pow(u[2],2)*u[1])/2 + (pow(u[1],2)*u[3])/2)));
		polymaster[5] += ((df(x,4)*pow(u[0],x-4)*pow(u[1],3)*u[2])/6) + (df(x,5)*pow(u[0],x-5)*pow(u[1],5)/120);
		return polymaster[5];
	}	
	else if(order == 6)
	{
		polymaster[6] = x*pow(u[0],x-1)*u[6] + (x*(x-1)*pow(u[0],x-2)*((u[3]*u[3]/2)+u[2]*u[4]+u[1]*u[5])) + df(x,3)*pow(u[0],x-3)*(pow(u[2],3)/6 + u[1]*u[2]*u[3] + u[1]*u[1]*u[4]/2);
		polymaster[6] += df(x,4)*pow(u[0],x-4)*((u[1]*u[1]*u[2]*u[2]/4) + (pow(u[1],3)*u[3]/6)) + (df(x,5)*pow(u[0],x-5)*pow(u[1],4)*u[2]/24) + (df(x,6)*pow(u[0],x-6)*pow(u[1],6)/720);
		return polymaster[6];
	}
	else if(order == 7)
	{
		polymaster[7] = x*pow(u[0],x-1)*u[7] + (x*(x-1)*pow(u[0],x-2)*((u[3]*u[4])+(u[2]*u[5])+(u[1]*u[6]))) + (df(x,3)*pow(u[0],x-3)*((u[2]*u[2]*u[3]/2)+(u[1]*u[3]*u[3]/2)+(u[1]*u[1]*u[5]/2)+(u[1]*u[2]*u[4])));
		polymaster[7] += df(x,4)*pow(u[0],x-4)*((u[1]*u[1]*u[2]*u[3]/2)+(pow(u[2],3)*u[1]/6)+(pow(u[1],3)*u[4]/6)) + (df(x,5)*pow(u[0],x-5)*((pow(u[1],4)*u[3]/24)+(pow(u[1],3)*u[2]*u[2]/12)));
		polymaster[7] += df(x,6)*pow(u[0],x-6)*pow(u[1],5)*u[2]/120 + df(x,7)*pow(u[0],x-7)*pow(u[1],7)/5040;
		return polymaster[7];
	}
	else if(order == 8)
	{
		polymaster[8] = x*pow(u[0],x-1)*u[8] + x*(x-1)*pow(u[0],x-2)*((u[4]*u[4]/2)+(u[3]*u[5])+(u[2]*u[6])+(u[1]*u[7])) + df(x,3)*pow(u[0],x-3)*((u[3]*u[3]*u[2]/2)+(u[2]*u[2]*u[4]/2)+(u[1]*u[1]*u[6]/2)+(u[1]*u[3]*u[4])+(u[1]*u[2]*u[5])); 
		polymaster[8] += df(x,4)*pow(u[0],x-4)*((u[1]*u[2]*u[2]*u[3]/2)+(pow(u[1],3)*u[5]/6)+(pow(u[2],4)/24)+(u[1]*u[1]*u[3]*u[3]/4)+(u[1]*u[1]*u[2]*u[4]/2)) + df(x,5)*pow(u[0],x-5)*((pow(u[2],3)*pow(u[1],2)/12)+(pow(u[1],4)*u[4]/24)+(pow(u[1],3)*u[2]*u[3]/6));
		polymaster[8] += df(x,6)*pow(u[0],x-6)*((pow(u[1],4)*u[2]*u[2]/48)+(pow(u[1],5)*u[3]/120)) + df(x,7)*pow(u[0],x-7)*pow(u[1],6)*u[2]/720 + df(x,8)*pow(u[0],x-8)*pow(u[1],8)/40320;
		return polymaster[8];
	}
	else if(order == 9)
	{
		polymaster[9] = x*pow(u[0],x-1)*u[9] + x*(x-1)*pow(u[0],x-2)*((u[4]*u[5])+(u[3]*u[6])+(u[2]*u[7])+(u[1]*u[8])) + df(x,3)*pow(u[0],x-3)*((pow(u[3],3)/6)+(u[2]*u[2]*u[5]/2)+(u[1]*u[4]*u[4]/2)+(u[7]*u[1]*u[1]/2)+(u[2]*u[3]*u[4])+(u[1]*u[3]*u[5])+(u[1]*u[2]*u[6]));
		polymaster[9] += df(x,4)*pow(u[0],x-4)*((u[1]*u[2]*u[3]*u[3]/2)+(pow(u[2],3)*u[3]/6)+(pow(u[1],3)*u[6]/6)+(u[1]*u[2]*u[2]*u[4]/2)+(u[1]*u[1]*u[3]*u[4]/2)+(u[1]*u[1]*u[2]*u[5]/2));
		polymaster[9] += df(x,5)*pow(u[0],x-5)*((pow(u[2],4)*u[1]/24)+(pow(u[1],4)*u[5]/24)+(pow(u[1],3)*u[3]*u[3]/12)+(u[1]*u[1]*u[2]*u[2]*u[3]/4)+(pow(u[1],3)*u[2]*u[4]/6)) + df(x,6)*(pow(u[0],x-6))*((pow(u[1],3)*pow(u[2],3)/36)+(pow(u[1],4)*u[2]*u[3]/24)+(pow(u[1],5)*u[4]/120)); 
		polymaster[9] += df(x,7)*pow(u[0],x-7)*((pow(u[1],5)*u[2]*u[2]/240)+(pow(u[1],6)*u[3]/720)) + df(x,8)*pow(u[0],x-8)*pow(u[1],7)*u[2]/5040 + df(x,9)*pow(u[1],9)/362880;
		return polymaster[9];
	}
}

double ADM(double temp, double pres, double yg, double yn, double yco, double yhs, double tpc, double ppc, double tpr, double ppr)
{
	double tr;
	if(yg!=-1)
	{
		tpc = 326 + 315.7*(yg-0.5) - (240*yn) - (83.3*yco) + (133.3*yhs);
		ppc = 678 - 50*(yg-0.5) - (206.7*yn) + (440*yco) + (606.7*yhs);
		tpc*=0.55555555556;//Since tpc is being calculated in rankine, we are converting it to kelvin(K).
		ppc*=0.006894757293;//Since ppc is being calculated in psia, we are converting it to MegaPascal(MP).
		tr = (tpc/temp);	
		ppr = pres/ppc;
	}
	else if(tpc!=-1)
	{
		tr = (tpc/temp);	
		ppr = pres/ppc;
	}
	else if(tpr!=-1)
	{
		tr=(1/tpr);
	}
	
	double A,B,C,D;
	A = (0.06125*tr)*(pow(e,(-1.2*(1-tr)*(1-tr))));
	B = tr*(14.76-(9.76*tr)+(4.58*tr*tr));
	C = tr*(90.7-(242.2*tr)+(42.4*tr*tr));
	D = 2.18 + (2.82*tr);
	
	double adompoly[10][8],term[8],first[8],y[11]={0,0,0,0,0,0,0,0,0,0,0},denom=(3*A*ppr+1),exp[8]={5,4,3,2,D+3,D+2,D+1,D},yf=0;
	y[0] = A*ppr/denom;
	first[0]=-B/denom;
	first[1]=(3*B+1)/denom;
	first[2]=(-(A*ppr)-(3*B)-1)/denom;
	first[3]=((3*A*ppr)+B-1)/denom;
	first[4]=C/denom;
	first[5]=-3*C/denom;
	first[6]=3*C/denom;
	first[7]=-C/denom;
	
	for(int i=0;i<10;i++)
	{
		for(int j=0;j<8;j++)
		{
			adompoly[i][j] = masterpoly(exp[j],y,i);
		}
		
		for(int j=0;j<8;j++)
			term[j] = first[j]*adompoly[i][j];
		for(int j=0;j<8;j++)
			y[i+1]+=term[j];
	}
	
	for(int i=0;i<11;i++)
	{
		yf+=y[i];
	}
	cout<<"Value of Y from ADM: "<<yf<<"\n";	
	cout<<"Value of z from ADM: "<<A*ppr/yf<<"\n\n";
	return (A*ppr)/yf;
}

double ADM_Shanks(double temp, double pres, double yg, double yn, double yco, double yhs, double tpc, double ppc, double tpr, double ppr)
{
	double tr;
	if(yg!=-1)
	{
		tpc = 326 + 315.7*(yg-0.5) - (240*yn) - (83.3*yco) + (133.3*yhs);
		ppc = 678 - 50*(yg-0.5) - (206.7*yn) + (440*yco) + (606.7*yhs);
		tpc*=0.55555555556;//Since tpc is being calculated in rankine, we are converting it to kelvin(K).
		ppc*=0.006894757293;//Since ppc is being calculated in psia, we are converting it to MegaPascal(MP).
		tr = (tpc/temp);	
		ppr = pres/ppc;
	}
	else if(tpc!=-1)
	{
		tr = (tpc/temp);	
		ppr = pres/ppc;
	}
	else if(tpr!=-1)
	{
		tr=(1/tpr);
	}
	
	double A,B,C,D;
	A = (0.06125*tr)*(pow(e,(-1.2*(1-tr)*(1-tr))));
	B = tr*(14.76-(9.76*tr)+(4.58*tr*tr));
	C = tr*(90.7-(242.2*tr)+(42.4*tr*tr));
	D = 2.18 + (2.82*tr);
	
	double adompoly[10][8],term[8],first[8],y[11]={0,0,0,0,0,0,0,0,0,0,0},denom=(3*A*ppr+1),exp[8]={5,4,3,2,D+3,D+2,D+1,D};
	y[0] = A*ppr/denom;
	first[0]=-B/denom;
	first[1]=(3*B+1)/denom;
	first[2]=(-(A*ppr)-(3*B)-1)/denom;
	first[3]=((3*A*ppr)+B-1)/denom;
	first[4]=C/denom;
	first[5]=-3*C/denom;
	first[6]=3*C/denom;
	first[7]=-C/denom;
	
	for(int i=0;i<5;i++)
	{
		for(int j=0;j<8;j++)
		{
			adompoly[i][j] = masterpoly(exp[j],y,i);
		}
		
		for(int j=0;j<8;j++)
			term[j] = first[j]*adompoly[i][j];
		for(int j=0;j<8;j++)
			y[i+1]+=term[j];
	}
	
	double ysum[5]={0,0,0,0,0},yshank[4];
	for(int i=0;i<5;i++)
	{
		for(int j=0;j<i+1;j++)
		{
			ysum[i]+=y[j];
		}
	}
	for(int i=0;i<3;i++)
	{
		yshank[i] = (ysum[i]*ysum[i+2] - ysum[i+1]*ysum[i+1])/(ysum[i+2]-2*ysum[i+1]+ysum[i]);
	}
	yshank[3] = (yshank[0]*yshank[2] - yshank[1]*yshank[1])/(yshank[2]-2*yshank[1]+yshank[0]);	
	cout<<"Value of Y from ADM-Shanks method: "<<yshank[3]<<"\n";
	cout<<"Value of z from ADM-Shanks method: "<<A*ppr/yshank[3]<<"\n\n";
	return (A*ppr)/yshank[3];
}

double Newton(double temp, double pres, double yg, double yn, double yco, double yhs, double tpc, double ppc, double tpr, double ppr,double ig)
{
	double tr;
	if(yg!=-1)
	{
		tpc = 326 + 315.7*(yg-0.5) - (240*yn) - (83.3*yco) + (133.3*yhs);
		ppc = 678 - 50*(yg-0.5) - (206.7*yn) + (440*yco) + (606.7*yhs);
		tpc*=0.55555555556;//Since tpc is being calculated in rankine, we are converting it to kelvin(K).
		ppc*=0.006894757293;//Since ppc is being calculated in psia, we are converting it to MegaPascal(MP).
		tr = (tpc/temp);	
		ppr = pres/ppc;
	}
	else if(tpc!=-1)
	{
		tr = (tpc/temp);	
		ppr = pres/ppc;
	}
	else if(tpr!=-1)
	{
		tr=(1/tpr);
	}	
	
	double A,B,C,D;
	A = (0.06125*tr)*(pow(e,(-1.2*(1-tr)*(1-tr))));
	B = tr*(14.76-(9.76*tr)+(4.58*tr*tr));
	C = tr*(90.7-(242.2*tr)+(42.4*tr*tr));
	D = 2.18 + (2.82*tr);
	
	double error,xprev,xnew,fx,fxi,counter=0,cons;

	error=100;
	
	xnew = ig;
	xprev = xnew;
	counter=1;
	while(counter<=10)
	{
		fx = ((xnew*( 1 + xnew + pow(xnew,2) - pow(xnew,3)))/(pow(1-xnew,3))) - (A*ppr) - (B*xnew*xnew) + (C*pow(xnew,D));
		fxi = ((pow(xnew,4)-(4*pow(xnew,3))+(4*xnew*xnew)+(4*xnew)+1)/pow(xnew-1,4)) - (2*B*xnew) + (D*C*pow(xnew,D-1));
		xnew -= (fx/fxi);
		xprev = xnew;
		counter+=1;
		cons=xnew;
	}
	cout<<"\nInital guess: "<<ig<<endl;
	cout<<"Value of Y from Newton Raphson method at end of 10th iteration: "<<cons<<"\n";
	cout<<"Value of z from Newton Raphson method at end of 10th iteration: "<<A*ppr/cons<<"\n\n";

	return (A*ppr)/cons;
}

int main()
{
	double temper[3]={337.872222,355.372,310.928},press[3]={13.789514,34.473786,6.894757},yg[2]={0.7,0.65},yn[2]={0.05,0.1},yc[2]={0.05,0.08},yh[2]={0.02,0.02};
	double tr[4]={1.3,1.5,1.7,2},pr[7]={0.5,1.5,2.5,3.5,4.5,5.5,6.5},count=4,a,zn,za,zs;
	
	auto start = chrono::high_resolution_clock::now();
	auto st = chrono::high_resolution_clock::now();
	a=Newton(temper[0],press[0],yg[0],yn[0],yc[0],yh[0],-1,-1,-1,-1,0.8);
	auto stop = chrono::high_resolution_clock::now();
	double t1 =  chrono::duration_cast<chrono::nanoseconds>(stop - st).count();
	
	a=ADM(temper[0],press[0],yg[0],yn[0],yc[0],yh[0],-1,-1,-1,-1);
	auto stopa = chrono::high_resolution_clock::now();
	double t2 =  chrono::duration_cast<chrono::nanoseconds>(stopa - stop).count();
	
	a=ADM_Shanks(temper[0],press[0],yg[0],yn[0],yc[0],yh[0],-1,-1,-1,-1);
	auto stopb = chrono::high_resolution_clock::now();
	double t3 =  chrono::duration_cast<chrono::nanoseconds>(stopb - stopa).count();
	
	auto end = chrono::high_resolution_clock::now();
	t1*=pow(10,-9);
	t2*=pow(10,-9);
	t3*=pow(10,-9);
	
	cout<<"For test case 1:\n";
	cout<<"Time taken in Newton-Raphson is :"<<t1<<"\n";
	cout<<"Time taken in ADM is :"<<t2<<"\n";
	cout<<"Time taken in ADM-Shanks is :"<<t3<<"\n\n";	
}
