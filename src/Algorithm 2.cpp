
#include<memory.h>
#include<vector>
#include<set>
#include<map>
#include<algorithm>
#include<iostream>
#include<fstream>
#include<time.h>
#include<ilcplex/cplex.h>

#define INFTY 10000
#define EPSILON 1E-5

using namespace std;

int m, n, cnt;

short* c = NULL;
short* maxrhs = NULL;
short* maxx = NULL;
double* indinc = NULL; // index increment vector

int nnegative; // number of columns with negative elements
int** matrix_s = NULL;

//////////////////////////////////////////////////////////////////////////

class point
{
public:
	short* vec;
	short valfun;
	double uniqueID;
	bool minimal; // equals 1 if point is lsm and <= maxrhs

	point() { vec = NULL; valfun = 0; uniqueID = 0.0; minimal = false; }

	~point() { if (vec != NULL) { delete[] vec; vec = NULL; } }

};

bool sortPointbyID(point* i, point* j) { return (i->uniqueID < j->uniqueID); }

vector<point> A;
vector<point*> lsm[2];

double src(const point* p) 
{
	double s = 0;
	for (int i = 0; i < m; i++) { s += p->vec[i] * indinc[i]; }
	return s;
}

short value(point* p, int cntt)
{
	// Compute z_k(p) over \bar{B}_k

	int i;
	bool smaller;
	bool stop, feasible;
	short vf = -INFTY;

	stop = false;

	for (vector<point*>::iterator it = lsm[cntt].begin(); it != lsm[cntt].end(); it++)
	{
		if ((*it)->valfun <= vf) { continue; }

		smaller = true;
		feasible = true;

		for (i = 0; i < m; i++)
		{
			if (p->vec[i] > (*it)->vec[i])
			{
				smaller = false;
			}

			if (p->vec[i] < (*it)->vec[i])
			{
				if (smaller == true)
				{
					stop = true;
				}

				feasible = false;

				break;
			}
		}

		if (stop == true) { break; }

		if (feasible == true) { vf = (*it)->valfun; }
	}

	return vf;
}

vector<point*>::iterator FindPoint(point* p, int cntt, bool* found)
{
	int i;
	bool psmaller;
	vector<point*>::iterator It;
	pair<vector<point*>::iterator, vector<point*>::iterator> bookends;

	*found = false;

	bookends = equal_range(lsm[cntt].begin(), lsm[cntt].end(), p, sortPointbyID);

	for (It = bookends.first; It != bookends.second; It++)
	{
		psmaller = false;

		for (i = 0; i < m; i++)
		{
			if (p->vec[i] < (*It)->vec[i])
			{
				psmaller = true;
				break;
			}

			if (p->vec[i] > (*It)->vec[i])
			{
				break;
			}
		}

		if (psmaller == true)
		{
			break;
		}

		if (i == m)
		{
			*found = true;
			break;
		}
	}

	return It;
}

void StorePoint(point* p, int cntt)
{
	int i;
	bool psmaller;
	vector<point*>::iterator It;
	pair<vector<point*>::iterator, vector<point*>::iterator> bookends;

	bookends = equal_range(lsm[cntt].begin(), lsm[cntt].end(), p, sortPointbyID);

	for (It = bookends.first; It != bookends.second; It++)
	{
		psmaller = false;

		for (i = 0; i < m; i++)
		{
			if (p->vec[i] < (*It)->vec[i])
			{
				psmaller = true;
				lsm[cntt].insert(It, p);
				break;
			}

			if (p->vec[i] > (*It)->vec[i])
			{
				break;
			}
		}

		if (psmaller == true || i == m)
		{
			break;
		}
	}

	if (It == bookends.second)
	{
		lsm[cntt].insert(It, p);
	}
}

void update_positive(int col)
{
	int i, j, t;
	int u_k;
	short nowvalue; // z_{k}(\beta)
	short prevalue;
	bool feasible;
	bool existence, existence2;
	vector<point*>::iterator pos, pos2, it;

	t = 0;
	point* zero = new point;
	zero->vec = new short[m];
	for (i = 0; i < m; i++) { zero->vec[i] = 0; }
	zero->valfun = 0;
	zero->uniqueID = 0;
	zero->minimal = true;
	StorePoint(zero, cnt);

	t = 1;
	point* a_k = new point;
	a_k->vec = new short[m];
	for (i = 0; i < m; i++) { a_k->vec[i] = A[col].vec[i]; }
	a_k->valfun = c[col];
	a_k->uniqueID = A[col].uniqueID;

	feasible = true;
	for (i = 0; i < m; i++) { if (a_k->vec[i] > maxrhs[i]) { feasible = false; break; } }
	if (feasible == true) { a_k->minimal = true; }

	StorePoint(a_k, cnt);

	t = 2;

	while (t <= maxx[col])
	{
		point* ta_k = new point;
		ta_k->vec = new short[m];
		for (i = 0; i < m; i++) { ta_k->vec[i] = t * A[col].vec[i]; }

		if (value(ta_k, 1 - cnt) >= t * c[col]) { break; }

		ta_k->valfun = t * c[col];
		ta_k->uniqueID = t * A[col].uniqueID;

		feasible = true;
		for (i = 0; i < m; i++) { if (ta_k->vec[i] > maxrhs[i]) { feasible = false; break; } }
		if (feasible == true) { ta_k->minimal = true; }

		StorePoint(ta_k, cnt);

		t++;
	}

	u_k = t - 1;

	it = lsm[1 - cnt].begin();

	while (it != lsm[1 - cnt].end())
	{
		if ((*it)->uniqueID == 0)
		{
			it++;
			continue;
		}

		// check the old ones

		point* now = new point;
		now->vec = new short[m];
		for (i = 0; i < m; i++) { now->vec[i] = (*it)->vec[i]; }
		now->uniqueID = (*it)->uniqueID;
		now->valfun = (*it)->valfun;

		feasible = true;
		for (j = 0; j < m; j++) { if (now->vec[j] > maxrhs[j]) { feasible = false; break; } }

		if (feasible == false)
		{
			it++;
			continue;
		}

		point* pre = new point;
		pre->vec = new short[m];
		// pre = now - A[col];
		for (i = 0; i < m; i++) { pre->vec[i] = now->vec[i] - A[col].vec[i]; }
		pre->uniqueID = now->uniqueID - A[col].uniqueID;

		pos = FindPoint(now, cnt, &existence);

		if (existence == true && (*pos)->minimal == true)
		{
			it++;
			continue;
		}

		pre->valfun = value(pre, cnt); // z_{k}(\beta - a_k)

		nowvalue = max(pre->valfun + c[col], (int)(*it)->valfun);

		now->valfun = nowvalue;

		if ((*it)->valfun > pre->valfun + c[col])
		{
			if (existence == true)
			{
				(*pos)->valfun = now->valfun;
				(*pos)->minimal = true;
			}
			else
			{
				now->minimal = true;
				lsm[cnt].insert(pos, now);
			}
		}
		else
		{
			pos2 = FindPoint(pre, cnt, &existence2);

			if (existence2 == true)
			{
				if (existence == true)
				{
					(*pos)->valfun = now->valfun;
					(*pos)->minimal = true;
				}
				else
				{
					now->minimal = true;
					lsm[cnt].insert(pos, now);
				}
			}
			else
			{
				if (existence == true)
				{
					lsm[cnt].erase(pos);
				}

				it++;
				continue;
			}
		}

		// find new ones

		for (i = 0; i < u_k; i++)
		{
			// larger = \beta + t a_k, nowvalue = z_{k}(\beta + (t-1)a_k)
			point* larger = new point;
			larger->vec = new short[m];
			for (j = 0; j < m; j++) { larger->vec[j] = now->vec[j] + (i + 1) * A[col].vec[j]; }
			larger->uniqueID = now->uniqueID + (i + 1) * (A[col].uniqueID);

			feasible = true;
			for (j = 0; j < m; j++) { if (larger->vec[j] > maxrhs[j]) { feasible = false; break; } }

			if (feasible == false)
			{
				break;
			}

			pos = FindPoint(larger, cnt, &existence);

			prevalue = value(larger, 1 - cnt); // z_{k-1}(\beta + t a_k)

			larger->valfun = max(nowvalue + c[col], (int)prevalue); // z_k(\beta + t a_k)

			if (prevalue < nowvalue + c[col])
			{
				if (existence == true)
				{
					(*pos)->valfun = larger->valfun;
					(*pos)->minimal = true;
				}
				else
				{
					larger->minimal = true;
					lsm[cnt].insert(pos, larger);
				}

				nowvalue = larger->valfun;
			}
			else
			{
				pos2 = FindPoint(larger, 1 - cnt, &existence2);

				if (existence2 == true)
				{
					if (existence == true)
					{
						(*pos)->valfun = larger->valfun;
						(*pos)->minimal = true;
					}
					else
					{
						larger->minimal = true;
						lsm[cnt].insert(pos, larger);
					}

					nowvalue = larger->valfun;
				}
				else
				{
					if (existence == true)
					{
						lsm[cnt].erase(pos);
					}
					
					break;
				}
			}
		}

		it++;
	}
}

void update_negative(int col)
{
	int i, j, t;
	int u_k;
	bool feasible;
	bool existence, existence2;
	vector<point*>::iterator pos, pos2, it;

	t = 0;
	point* zero = new point;
	zero->vec = new short[m];
	for (i = 0; i < m; i++) { zero->vec[i] = 0; }
	zero->valfun = 0;
	zero->uniqueID = 0;
	zero->minimal = true;
	StorePoint(zero, cnt);

	t = 1;
	point* a_k = new point;
	a_k->vec = new short[m];
	for (i = 0; i < m; i++) { a_k->vec[i] = A[col].vec[i]; }
	a_k->valfun = c[col];
	a_k->uniqueID = A[col].uniqueID;

	feasible = true;
	for (i = 0; i < m; i++) { if (a_k->vec[i] > maxrhs[i]) { feasible = false; break; } }
	if (feasible == true) { a_k->minimal = true; }

	StorePoint(a_k, cnt);

	t = 2;

	while (t <= maxx[col])
	{
		point* ta_k = new point;
		ta_k->vec = new short[m];
		for (i = 0; i < m; i++) { ta_k->vec[i] = t * A[col].vec[i]; }

		if (value(ta_k, 1 - cnt) >= t * c[col]) { break; }

		ta_k->valfun = t * c[col];
		ta_k->uniqueID = t * A[col].uniqueID;

		feasible = true;
		for (i = 0; i < m; i++) { if (ta_k->vec[i] > maxrhs[i]) { feasible = false; break; } }
		if (feasible == true) { ta_k->minimal = true; }

		StorePoint(ta_k, cnt);

		t++;
	}

	u_k = t - 1;

	it = lsm[1 - cnt].begin();

	while (it != lsm[1 - cnt].end())
	{
		if ((*it)->uniqueID == 0)
		{
			it++;
			continue;
		}

		// check the old ones

		point* now = new point;
		now->vec = new short[m];
		for (i = 0; i < m; i++) { now->vec[i] = (*it)->vec[i]; }
		now->uniqueID = (*it)->uniqueID;
		now->valfun = (*it)->valfun;

		pos = FindPoint(now, cnt, &existence);

		if (existence == true && (*pos)->minimal == true)
		{
			it++;
			continue;
		}

		if (existence == true)
		{
			if ((*pos)->valfun < now->valfun)
			{
				(*pos)->valfun = now->valfun;
			}
			else
			{
				it++;
				continue;
			}
		}
		else
		{
			lsm[cnt].insert(pos, now);
		}

		// check the new ones

		for (i = 0; i < u_k; i++)
		{
			// larger = \beta + t a_k, nowvalue = z_{k}(\beta + (t-1)a_k)
			point* larger = new point;
			larger->vec = new short[m];
			for (j = 0; j < m; j++) { larger->vec[j] = now->vec[j] + (i + 1) * A[col].vec[j]; }
			larger->uniqueID = now->uniqueID + (i + 1) * (A[col].uniqueID);
			larger->valfun = now->valfun + (i + 1) * c[col];

			pos2 = FindPoint(larger, cnt, &existence2);

			if (existence2 == true)
			{
				if ((*pos2)->valfun < larger->valfun)
				{
					(*pos2)->valfun = larger->valfun;
				}
				else
				{
					break;
				}
			}
			else
			{
				lsm[cnt].insert(pos2, larger);
			}
		}

		it++;
	}
}

void readdata(string datafile_name)
{
	int i, j;

	ifstream datafile;
	datafile.open(datafile_name, ifstream::in);
	if (!datafile.is_open()) {
		cout << "Problem: data file is not available!" << endl;
		exit(1);
	}

	// m: number of rows
	// n: number of columns

	datafile >> m >> n;

	maxrhs = new short[m];
	for (i = 0; i < m; i++) { datafile >> maxrhs[i]; }

	//====================================================
	for (i = 0; i < m; i++) { maxrhs[i] = maxrhs[i] + 0; }
	//====================================================

	c = new short[n];
	for (j = 0; j < n; j++) { datafile >> c[j]; }

	A.resize(n);

	for (j = 0; j < n; j++) { A[j].vec = new short[m]; }

	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++) { datafile >> A[j].vec[i]; }
	}

	datafile.close();

	/*for (j = 0; j < n; j++) cout << c[j] << " ";
	cout << endl << endl;

	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++)
			cout << A[j].vec[i] << " ";
		cout << endl;
	}
	cout << endl;*/

}

void output()
{
	int k = -1;

	ofstream resultfile;

	resultfile.open("out_cbc.txt", ofstream::out);

	resultfile << "size = " << lsm[cnt].size() << endl << endl;

	for (vector<point*>::iterator it = lsm[cnt].begin(); it != lsm[cnt].end(); it++)
	{
		k++;

		resultfile << "k = " << k << ": ";

		for (int i = 0; i < m; i++) { resultfile << (*it)->vec[i] << " "; }

		resultfile << "  v = " << (*it)->valfun << endl;
	}

	resultfile.close();
}

int createmasterproblem(CPXENVptr env, CPXLPptr lp)
{
	//  max c'x
	// s.t. Ax <= beta, x>=0 integer

	int i, j;
	int status = 0;

	double* obj = new double[n];
	double* lb = new double[n];
	double* ub = new double[n];
	char* ctype = new char[n];

	int rmatbeg;
	int* rmatind = new int[n];
	double* rmatval = new double[n];
	double rhs;
	char sense;

	/* Change optimization sense */
	status = CPXchgobjsen(env, lp, CPX_MAX);
	if (status) {
		fprintf(stderr, "Error in CPXchgobjsen in createmasterproblem(), status = %d.\n", status);
		goto TERMINATE;
	}

	/* Change problem type */
	status = CPXchgprobtype(env, lp, CPXPROB_MILP);
	if (status) {
		fprintf(stderr, "Error in CPXchgprobtype in createmasterproblem(), status = %d.\n", status);
		goto TERMINATE;
	}

	/* create variables and add columns */
	for (i = 0; i < n; i++) {
		obj[i] = c[i];
		lb[i] = 0;
		ub[i] = CPX_INFBOUND;
		ctype[i] = 'I';
	}
	status = CPXnewcols(env, lp, n, obj, lb, ub, ctype, NULL);
	if (status) {
		fprintf(stderr, "Error in CPXnewcols in createmasterproblem(), status = %d.\n", status);
		goto TERMINATE;
	}

	/* add constraints */
	rmatbeg = 0;
	sense = 'L';
	for (i = 0; i < n; i++) { rmatind[i] = i; }
	for (i = 0; i < m; i++) {
		rhs = maxrhs[i];
		for (j = 0; j < n; j++) { rmatval[j] = A[j].vec[i]; }
		status = CPXaddrows(env, lp, 0, 1, n, &rhs, &sense, &rmatbeg, rmatind, rmatval, NULL, NULL);
		if (status) {
			fprintf(stderr, "Error in CPXaddrows in createmasterproblem(), status = %d.\n", status);
			goto TERMINATE;
		}
	}

	//CPXwriteprob(env, lp, "masterproblem.lp", NULL);

TERMINATE:

	delete[] obj; obj = NULL;
	delete[] lb; lb = NULL;
	delete[] ub; ub = NULL;
	delete[] ctype; ctype = NULL;

	delete[] rmatind; rmatind = NULL;
	delete[] rmatval; rmatval = NULL;

	return status;
}

int FilterPoint()
{
	// Filter those dominated monoid elements

	int i;
	int counter = 0;
	bool agb, bga, remark;
	vector<point*>::iterator It, It2, tempIt;

	for (It = lsm[cnt].begin(); It != lsm[cnt].end();)
	{
		// if ((*It)->minimal == true) { It++; continue; } 
		// What if It dominates some other vector? Those dominated vectors are not smaller than It. 
		// This can be wrong. For example, beta <= max rhs, beta-a_k not <= max rhs, but beta+a_k <= max rhs.
		// In this case, beta+a_k may be dominated by a lsm vector.

		tempIt = It;
		tempIt++;
		remark = false;

		for (It2 = tempIt; It2 != lsm[cnt].end();)
		{
			//if ((*It2)->minimal == true) { It2++; continue; } // What if It2 dominates some other vector?

			agb = false;
			bga = false;

			for (i = 0; i < m; i++)
			{
				if ((agb == false) || (bga == false)) {
					if ((*It)->vec[i] > (*It2)->vec[i])
						agb = true;
					else if ((*It)->vec[i] < (*It2)->vec[i])
						bga = true;
				}
				else { break; }
			}

			if ((agb == true) && (bga == false))
			{
				if ((*It)->valfun <= (*It2)->valfun)
				{
					It = lsm[cnt].erase(It);
					counter++;
					remark = true;
					break;
				}
				else { It2++; }
			}
			else if ((agb == false) && (bga == true))
			{
				if ((*It)->valfun >= (*It2)->valfun)
				{
					It2 = lsm[cnt].erase(It2);
					counter++;
				}
				else { It2++; }
			}
			else { It2++; }
		}

		if (remark == false) { It++; }
	}

	return (counter);
}

void FilterLargePoint(int col)
{
	// Filter those monoid elements that are out of bound and not possibly become smaller afterwards

	int i;
	bool feasible;
	vector<point*>::iterator it;

	for (it = lsm[cnt].begin(); it != lsm[cnt].end(); )
	{
		feasible = true;

		for (i = 0; i < m; i++)
		{
			if ((*it)->vec[i] + matrix_s[i][col] > maxrhs[i])
			{
				feasible = false;
				break;
			}
		}

		if (feasible == true)
		{
			it++;
		}
		else
		{
			/*for (i = 0; i < m; i++) { cout << (*it)->vec[i] << " "; }
			cout << endl;*/
			it = lsm[cnt].erase(it);
		}
	}
}


int main(int argc, char* argv[])
{
	int i, j, k;
	int status = 0;
	int deletedcolumns;
	int t;
	bool feasible;
	
	int solstat;
	double objval;
	double* solution;
	int* indices;
	double* newobj;
	char* lu;
	double* newbd;
	double* newrhs;

	int* minax; // minimum value that each dimension of any lsm can take
	short* nnegative_index = NULL;

	vector<point*>::iterator it;

	CPXENVptr env = NULL;
	CPXLPptr  lp = NULL;

	time_t tstart_global, tstart_local, optend;
	double time_limit; // runtime limit
	double time_preprocessing; // pre-processing time
	double time_negative; // runtime for processing nonpositive columns
	double time_negative_filter; // runtime for filtering non-lsm elements when processing nonpositive columns
	double time_filter; // runtime for filtering elements not <= maxrhs after processing nonpositive columns
	double time_positive; // runtime for processing nonnegative columns
	double time_columnbycolumn; // total runtime of column by column algorithm
	double time_elapse;

	tstart_global = clock(); // global time start

	string folder = "E:\\Data\\Andy instances\\Instances_negative\\";
	// string instance = "Instance-4.txt";
	string instance = argv[1];
	string instancedir = folder + instance;
	string resultdir = folder + "Result_cbc.txt";
	string exceldir = folder + "Result_cbc_excel.txt";

	ofstream resultfile, excelfile;
	resultfile.open(resultdir, ofstream::out | ofstream::app);
	excelfile.open(exceldir, ofstream::out | ofstream::app);

	// time_limit = 21600;
	time_limit = atof(argv[2]);

	cout << endl << "Problem: " << instance << endl << endl;
	resultfile << "Problem: " << instance << endl << endl;
	
	readdata(instancedir);

	solution = new double[n];
	maxx = new short[n];
	indices = new int[n];
	newobj = new double[n];
	lu = new char[n];
	newbd = new double[n];
	newrhs = new double[m];
	minax = new int[m];


	//////////////////////////////////////////////////////////////////////////
	// Check the number of columns with negative elements

	nnegative = 0;

	for (j = 0; j < n; j++)
	{
		bool positive = true;

		for (i = 0; i < m; i++)
		{
			if (A[j].vec[i] < 0)
			{
				positive = false;
				nnegative++;
				break;
			}
		}

		if (positive == true)
		{
			break;
		}
	}


	//////////////////////////////////////////////////////////////////////////
	// Create master problem

	/* Init the CPLEX environment */
	env = CPXopenCPLEX(&status);
	if (env == NULL) {
		fprintf(stderr, "Failure in CPXopenCPLEX in main(), status = %d.\n", status);
		resultfile << "  Failure in CPXopenCPLEX in main(), status = " << status << endl;
		goto TERMINATE;
	}

	/* Create the master problem */
	lp = CPXcreateprob(env, &status, "masterproblem.lp");
	if (lp == NULL) {
		fprintf(stderr, "Failure in CPXcreateprob in main(), status = %d.\n", status);
		resultfile << "  Failure in CPXcreateprob in main(), status = " << status << endl;
		goto TERMINATE;
	}

	status = createmasterproblem(env, lp);
	if (status) {
		fprintf(stderr, "Failed to create the master problem in main().\n");
		resultfile << "  Failed to create the master problem in main()" << endl;
		goto TERMINATE;
	}

	for (j = 0; j < n; j++) { indices[j] = j; }


	//////////////////////////////////////////////////////////////////////////
	// Check if z(0) = 0

	/* Update the right-hand side vector of masterproblem */
	for (i = 0; i < m; i++) { newrhs[i] = 0; }

	status = CPXchgrhs(env, lp, m, indices, newrhs);
	if (status) {
		fprintf(stderr, "Error in CPXchgrhs in main(): status = %d\n", status);
		resultfile << "  Error in CPXchgrhs in main(): status = " << status << endl;
		goto TERMINATE;
	}

	/* Optimize the problem and obtain solution */
	status = CPXmipopt(env, lp);
	if (status) {
		fprintf(stderr, "Failed to optimize MIP in main().\n");
		resultfile << "  Failed to optimize MIP in main()" << endl;
		goto TERMINATE;
	}

	solstat = CPXgetstat(env, lp);

	if (solstat != 101)
	{
		if (solstat == 118)
		{
			fprintf(stderr, "z(0) != 0: CPXMIP_UNBOUNDED!  Exiting...\n");
			resultfile << "  z(0) != 0: CPXMIP_UNBOUNDED!" << endl;
		}
		else
		{
			printf("\nSolution status = %d\n", solstat);
			resultfile << "  Solution status = " << solstat << endl;
		}

		goto TERMINATE;
	}
	

	//////////////////////////////////////////////////////////////////////////
	// Find the maximum value that each variable can take with the max rhs

	/* Update the right-hand side vector of masterproblem */
	for (i = 0; i < m; i++) { newrhs[i] = maxrhs[i]; }

	status = CPXchgrhs(env, lp, m, indices, newrhs);
	if (status) {
		fprintf(stderr, "Error in CPXchgrhs in main(): status = %d\n", status);
		resultfile << "  Error in CPXchgrhs in main(): status = " << status << endl;
		goto TERMINATE;
	}
	
	for (j = 0; j < n; j++)
	{
		/* Update the objective function of masterproblem */
		for (k = 0; k < n; k++) { newobj[k] = 0; }
		newobj[j] = 1;

		status = CPXchgobj(env, lp, n, indices, newobj);
		if (status) {
			fprintf(stderr, "Error in CPXchgobj in main(): status = %d\n", status);
			resultfile << "  Error in CPXchgobj in main(): status = " << status << endl;
			goto TERMINATE;
		}

		/* Optimize the problem and obtain solution */
		status = CPXmipopt(env, lp);
		if (status) {
			fprintf(stderr, "Failed to optimize MIP in main().\n");
			resultfile << "  Failed to optimize MIP in main()." << endl;
			goto TERMINATE;
		}

		solstat = CPXgetstat(env, lp);

		if (solstat != 101)
		{
			if (solstat == 118)
			{
				fprintf(stderr, "x_k: CPXMIP_UNBOUNDED!  Exiting...\n");
				resultfile << "  x_k: CPXMIP_UNBOUNDED!" << endl;
			}
			else
			{
				printf("\nSolution status = %d\n", solstat);
				resultfile << "  Solution status = " << solstat << endl;
			}

			goto TERMINATE;
		}

		status = CPXgetobjval(env, lp, &objval);
		if (status) {
			fprintf(stderr, "No MIP objective value available.  Exiting...\n");
			resultfile << "  No MIP objective value available." << endl;
			goto TERMINATE;
		}

		maxx[j] = objval;

		/*printf("Solution value  = %f\n\n", objval);

		status = CPXgetx(env, lp, solution, 0, n - 1);
		if (status) {
			fprintf(stderr, "Failed to get optimal integer x.\n");
			goto TERMINATE;
		}

		for (j = 0; j < n; j++) { cout << solution[j] << " "; }
		cout << endl << endl;*/
	}

	cout << "  Max x: ";
	for (j = 0; j < n; j++) { cout << maxx[j] << " "; }
	cout << endl << endl;


	//////////////////////////////////////////////////////////////////////////
	// Find the minimum value that each dimension of any lsm can take with the max rhs

	/* Change optimization sense */
	status = CPXchgobjsen(env, lp, CPX_MIN);
	if (status) {
		fprintf(stderr, "Error in CPXchgobjsen in main(), status = %d.\n", status);
		resultfile << "  Error in CPXchgobjsen in main(), status = " << status << endl;
		goto TERMINATE;
	}

	for (i = 0; i < m; i++)
	{
		bool positive = true;

		for (j = 0; j < n; j++)
		{
			if (A[j].vec[i] < 0)
			{
				positive = false;
				break;
			}
		}

		if (positive == true) 
		{ 
			minax[i] = 0;
			continue; 
		}

		/* Update the objective function of masterproblem */
		for (j = 0; j < n; j++) { newobj[j] = A[j].vec[i]; }

		status = CPXchgobj(env, lp, n, indices, newobj);
		if (status) {
			fprintf(stderr, "Error in CPXchgobj in main(): status = %d\n", status);
			resultfile << "  Error in CPXchgobj in main(): status = " << status << endl;
			goto TERMINATE;
		}

		/* Optimize the problem and obtain solution */
		status = CPXmipopt(env, lp);
		if (status) {
			fprintf(stderr, "Failed to optimize MIP in main().\n");
			resultfile << "  Failed to optimize MIP in main()." << endl;
			goto TERMINATE;
		}

		status = CPXgetobjval(env, lp, &objval);
		if (status) {
			fprintf(stderr, "No MIP objective value available.  Exiting...\n");
			resultfile << "  No MIP objective value available." << endl;
			goto TERMINATE;
		}

		minax[i] = objval;
	}

	cout << "  Min a'x: ";
	for (i = 0; i < m; i++) { cout << minax[i] << " "; }
	cout << endl << endl;

	indinc = new double[m];
	indinc[m - 1] = 1;
	for (i = m - 2; i >= 0; i--) { indinc[i] = indinc[i + 1] * (maxrhs[i + 1] - minax[i + 1] + 1); }

	for (j = 0; j < n; j++) { A[j].uniqueID = src(&A[j]); }

	/* Restore the objective function of masterproblem */
	for (j = 0; j < n; j++) { newobj[j] = c[j]; }

	status = CPXchgobj(env, lp, n, indices, newobj);
	if (status) {
		fprintf(stderr, "Error in CPXchgobj in main(): status = %d\n", status);
		goto TERMINATE;
	}

	////////////////////////////////////////////////

	// Sort the non-positive columns in decreasing order of maxx[]

	nnegative_index = new short[nnegative];

	for (j = 0; j < nnegative; j++) { nnegative_index[j] = j; }

	for (j = 0; j < nnegative; j++)
	{
		for (i = 0; i < nnegative - 1 - j; i++)
		{
			if (maxx[i] < maxx[i + 1])
			{
				short t = maxx[i];
				maxx[i] = maxx[i + 1];
				maxx[i + 1] = t;

				t = nnegative_index[i];
				nnegative_index[i] = nnegative_index[i + 1];
				nnegative_index[i + 1] = t;

				t = c[i];
				c[i] = c[i + 1];
				c[i + 1] = t;

				point temp_point;
				temp_point.vec = new short[m];

				for (k = 0; k < m; k++) { temp_point.vec[k] = A[i].vec[k]; }
				temp_point.uniqueID = A[i].uniqueID;

				for (k = 0; k < m; k++) { A[i].vec[k] = A[i + 1].vec[k]; }
				A[i].uniqueID = A[i + 1].uniqueID;

				for (k = 0; k < m; k++) { A[i + 1].vec[k] = temp_point.vec[k]; }
				A[i + 1].uniqueID = temp_point.uniqueID;
			}
		}
	}

	matrix_s = new int* [m];
	for (i = 0; i < m; i++)
	{
		matrix_s[i] = new int[nnegative];
		for (j = 0; j < nnegative; j++) { matrix_s[i][j] = 0; }
	}

	for (j = nnegative - 2; j >= 0; j--)
	{
		for (i = 0; i < m; i++)
		{
			matrix_s[i][j] = matrix_s[i][j + 1];

			if (A[j + 1].vec[i] < 0)
			{
				matrix_s[i][j] += A[j + 1].vec[i] * maxx[j + 1];
			}
		}
	}

	/*for (i = 0; i < m; i++) {
		for (j = 0; j < nnegative; j++)
			cout << matrix_s[i][j] << " ";
		cout << endl;
	}
	cout << endl;*/

	/*for (j = 0; j < nnegative; j++) { cout << maxx[j] << ' '; }
	cout << endl;
	for (j = 0; j < nnegative; j++) { cout << nnegative_index[j] << ' '; }
	cout << endl;

	for (i = 0; i < m; i++) {
		for (j = 0; j < nnegative; j++)
			cout << A[j].vec[i] << " ";
		cout << endl;
	}
	cout << endl;

	for (j = 0; j < nnegative; j++) { cout << A[j].uniqueID << ' '; }
	cout << endl;*/

	////////////////////////////////////////////////

	optend = clock();
	time_preprocessing = (double)(optend - tstart_global) / CLOCKS_PER_SEC;


	//////////////////////////////////////////////////////////////////////////
	// Column by column begins
	
	// The first column

	tstart_local = clock(); // record the time for processing nonpositive columns

	cout << "  Number of nonpositive columns: " << nnegative << endl << endl;

	cout << "  Column: " << 0 << endl;

	deletedcolumns = 0;

	cnt = 0;

	t = 0;

	while (t <= maxx[0])
	{
		point* p = new point;

		p->vec = new short[m];
		for (i = 0; i < m; i++) { p->vec[i] = t * A[0].vec[i]; }

		p->valfun = t * c[0];
		p->uniqueID = t * A[0].uniqueID;

		feasible = true;
		for (i = 0; i < m; i++) { if (p->vec[i] > maxrhs[i]) { feasible = false; break; } }

		if (feasible == true) { p->minimal = true; }
		
		StorePoint(p, cnt);

		t++;
	}

	cout << "    Number of monoids generated: " << lsm[cnt].size() << endl;

	FilterLargePoint(0);

	cout << "    Number of monoids left: " << lsm[cnt].size() << endl;

	// The columns with negative elements

	for (k = 1; k < nnegative; k++)
	{
		cout << "  Column: " << k << endl;

		if (maxx[k] == 0 || value(&A[k], cnt) >= c[k])
		{
			cout << "    Column is eliminated!" << endl;
			deletedcolumns++;
			continue;
		}

		cnt = 1 - cnt;

		lsm[cnt].clear();

		update_negative(k);

		cout << "    Number of monoids generated: " << lsm[cnt].size() << endl;

		FilterLargePoint(k);

		cout << "    Number of monoids left: " << lsm[cnt].size() << endl;

		optend = clock();
		time_elapse = (double)(optend - tstart_global) / CLOCKS_PER_SEC;
		if (time_elapse > time_limit)
		{
			cout << endl << "  Runtime limit exceeded for nonpositive columns!" << endl << endl;
			resultfile << endl << "  Runtime limit exceeded for nonpositive columns!" << endl << endl;
			goto TERMINATE;
		}

		/*cout << endl << "column = " << k << endl << endl;
		for (set<point, cmp>::iterator it = lsm[cnt].begin(); it != lsm[cnt].end(); ++it)
		{
			for (i = 0; i < m; i++) { cout << it->vec[i] << " "; }

			cout << "  v = " << it->valfun << endl;
		}
		cout << endl;*/
	}

	optend = clock();
	time_negative = (double)(optend - tstart_local) / CLOCKS_PER_SEC;


	// Filter the monoid that is out of bound

	tstart_local = clock(); // record the time for filtering elements not <= maxrhs

	for (it = lsm[cnt].begin(); it != lsm[cnt].end(); )
	{
		if ((*it)->minimal == true) { it++; continue; }

		feasible = true;
		for (i = 0; i < m; i++) { if ((*it)->vec[i] > maxrhs[i]) { feasible = false; break; } }

		if (feasible == true) {
			it++;
		}
		else {
			it = lsm[cnt].erase(it);
		}
	}

	optend = clock();
	time_filter = (double)(optend - tstart_local) / CLOCKS_PER_SEC;

	cout << endl << "  Number of monoids left: " << lsm[cnt].size() << endl << endl;

	tstart_local = clock(); // record the time for filtering elements not lsm

	FilterPoint();

	optend = clock();
	time_negative_filter = (double)(optend - tstart_local) / CLOCKS_PER_SEC;

	cout << "  Number of lsm generated: " << lsm[cnt].size() << endl << endl;


	// The nonnegative columns

	tstart_local = clock(); // record the time for processing nonnegative columns

	for (k = nnegative; k < n; k++)
	{
		cout << "  Column: " << k << endl;

		if (maxx[k] == 0 || value(&A[k], cnt) >= c[k])
		{
			cout << "    Column is eliminated!" << endl;
			deletedcolumns++;
			continue;
		}

		cnt = 1 - cnt;

		lsm[cnt].clear();

		update_positive(k);

		cout << "    Number of lsm generated: " << lsm[cnt].size() << endl;

		optend = clock();
		time_elapse = (double)(optend - tstart_global) / CLOCKS_PER_SEC;
		if (time_elapse > time_limit)
		{
			cout << endl << "  Runtime limit exceeded for nonnegative columns!" << endl << endl;
			resultfile << endl << "  Runtime limit exceeded for nonnegative columns!" << endl << endl;
			goto TERMINATE;
		}

		/*cout << endl << "column = " << k << endl << endl;
		for (set<point, cmp>::iterator it = lsm[cnt].begin(); it != lsm[cnt].end(); ++it)
		{
			for (i = 0; i < m; i++) { cout << it->vec[i] << " "; }

			cout << "  v = " << it->valfun << endl;
		}
		cout << endl;*/
	}

	optend = clock();
	time_positive = (double)(optend - tstart_local) / CLOCKS_PER_SEC;
	

	// delete non-lsm vectors

	for (it = lsm[cnt].begin(); it != lsm[cnt].end(); )
	{
		if ((*it)->minimal == false)
		{
			it = lsm[cnt].erase(it);
		}
		else {
			it++;
		}
	}


	cout << endl << "  The number of columns: " << n << endl;
	cout << "  The number of columns deleted: " << deletedcolumns << endl;
	cout << "  The number of remaining columns: " << n - deletedcolumns << endl;
	cout << "  Number of lsm generated: " << lsm[cnt].size() << endl;

	resultfile << "  The number of columns: " << n << endl;
	resultfile << "  The number of columns deleted: " << deletedcolumns << endl;
	resultfile << "  The number of remaining columns: " << n - deletedcolumns << endl;
	resultfile << "  Number of lsm generated: " << lsm[cnt].size() << endl;

	excelfile << deletedcolumns << ' ' << lsm[cnt].size() << ' ';

	optend = clock();
	time_columnbycolumn = (double)(optend - tstart_global) / CLOCKS_PER_SEC;

	cout << "  Runtime for pre-processing: " << time_preprocessing << endl;
	cout << "  Runtime for processing nonpositive columns: " << time_negative << endl;
	cout << "  Runtime for filtering monoids out of bound: " << time_filter << endl;
	cout << "  Runtime for filtering non-lsm monoids: " << time_negative_filter << endl;
	cout << "  Runtime for processing nonnegative columns: " << time_positive << endl;
	cout << "  Total runtime for column by column: " << time_columnbycolumn << endl << endl;

	resultfile << "  Runtime for pre-processing: " << time_preprocessing << endl;
	resultfile << "  Runtime for processing nonpositive columns: " << time_negative << endl;
	resultfile << "  Runtime for filtering monoids out of bound: " << time_filter << endl;
	resultfile << "  Runtime for filtering non-lsm monoids: " << time_negative_filter << endl;
	resultfile << "  Runtime for processing nonnegative columns: " << time_positive << endl;
	resultfile << "  Total runtime for column by column: " << time_columnbycolumn << endl << endl;

	excelfile << time_preprocessing << ' ' << time_negative << ' ' << time_negative_filter << ' ';
	excelfile << time_filter << ' ' << time_positive << ' ' << time_columnbycolumn << ' ';

	excelfile << endl;

	//////////////////////////////////////////////////////////////////////////

TERMINATE:

	resultfile.close();
	excelfile.close();

	lsm[0].clear();
	lsm[1].clear();

	delete[] c; c = NULL;
	delete[] maxrhs; maxrhs = NULL;
	delete[] indinc; indinc = NULL;
	delete[] solution; solution = NULL;
	delete[] maxx; maxx = NULL;
	delete[] indices; indices = NULL;
	delete[] newobj; newobj = NULL;
	delete[] lu; lu = NULL;
	delete[] newbd; newbd = NULL;
	delete[] newrhs; newrhs = NULL; 
	delete[] minax; minax = NULL;
	delete[] nnegative_index; nnegative_index = NULL;

	for (i = 0; i < m; i++) { delete[] matrix_s[i]; matrix_s[i] = NULL; }
	delete[] matrix_s; matrix_s = NULL;

	if (lp != NULL) {
		int local_status = CPXfreeprob(env, &lp);
		if (local_status) {
			fprintf(stderr, "CPXfreeprob failed in main(), error code %d.\n", local_status);
			status = local_status;
		}
	}

	if (env != NULL) {
		int local_status = CPXcloseCPLEX(&env);
		if (local_status) {
			fprintf(stderr, "Could not close CPLEX environment in main(), status = %d.\n", local_status);
			status = local_status;
		}
	}

	//getchar();

	return 0;
}