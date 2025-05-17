
#include <memory.h>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <time.h>

#define INFTY 10000
#define EPSILON 1E-5

using namespace std;

int m, n, cnt = 0;

short* c = NULL;
short* maxrhs = NULL;
double* indinc = NULL; // index increment vector

double tlimit2; // running time limit for enumeration
time_t tstart;

//////////////////////////////////////////////////////////////////////////

class point
{
public:
	short* vec;
	short valfun;
	double uniqueID;

	point() { vec = NULL; valfun = 0; uniqueID = 0.0; }

	~point() { if (vec != NULL) { delete[] vec; vec = NULL; } }

};

bool sortPointbyID(point* i, point* j) { return (i->uniqueID < j->uniqueID); }

double src(const point *p) 
{
	double s = 0;
	for (int i = 0; i < m; i++) { s += p->vec[i] * indinc[i]; }
	return s;
}

point ub;
vector<point> A;
vector<point*> lsm[2];

short value(point* p, int cntt)
{
	int i;
	bool pequal;
	bool stop, feasible;
	short vf = -INFTY;

	stop = false;

	for (vector<point*>::iterator it = lsm[cntt].begin(); it != lsm[cntt].end(); it++)
	{
		if ((*it)->valfun <= vf) { continue; }

		pequal = true;
		feasible = true;

		for (i = 0; i < m; i++)
		{
			if (p->vec[i] > (*it)->vec[i])
			{
				pequal = false;
			}

			if (p->vec[i] < (*it)->vec[i])
			{
				if (pequal == true)
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

bool FindPoint(point* p, int cntt)
{
	int i;
	bool found = false;
	vector<point*>::iterator It;
	pair<vector<point*>::iterator, vector<point*>::iterator> bookends;

	bookends = equal_range(lsm[cntt].begin(), lsm[cntt].end(), p, sortPointbyID);

	for (It = bookends.first; It != bookends.second; It++)
	{
		for (i = 0; i < m; i++) {
			if ((*It)->vec[i] != p->vec[i]) { break; }
		}
		if (i == m) {
			found = true;
			break;
		}
	}

	return found;
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

void update(int col)
{
	int i, j, t;
	int u_k;
	short nowvalue;
	short prevalue;
	bool feasible;

	t = 0;
	point* zero = new point;
	zero->vec = new short[m];
	for (i = 0; i < m; i++) { zero->vec[i] = 0; }
	zero->valfun = 0;
	zero->uniqueID = 0;
	lsm[cnt].push_back(zero);

	t = 1;
	point* a_k = new point;
	a_k->vec = new short[m];
	for (i = 0; i < m; i++) { a_k->vec[i] = A[col].vec[i]; }
	a_k->valfun = c[col];
	a_k->uniqueID = A[col].uniqueID;
	lsm[cnt].push_back(a_k);

	t = 2;

	while (1)
	{
		point* ta_k = new point;
		ta_k->vec = new short[m];
		for (i = 0; i < m; i++) { ta_k->vec[i] = t * A[col].vec[i]; }

		feasible = true;
		for (i = 0; i < m; i++) { if (ta_k->vec[i] > ub.vec[i]) { feasible = false; break; } }

		if (feasible == false) { break; }

		if (value(ta_k, 1 - cnt) >= t * c[col]) { break; }

		ta_k->valfun = t * c[col];
		ta_k->uniqueID = t * A[col].uniqueID;
		lsm[cnt].push_back(ta_k);
		t++;
	}

	u_k = t - 1;

	vector<point*>::iterator it = lsm[1 - cnt].begin();

	while (it != lsm[1 - cnt].end())
	{
		if ((*it)->uniqueID == 0)
		{
			it++;
			continue;
		}

		// update the old ones
		point* now = new point;
		now->vec = new short[m];
		for (i = 0; i < m; i++) { now->vec[i] = (*it)->vec[i]; }
		now->uniqueID = (*it)->uniqueID;
		now->valfun = (*it)->valfun;

		point* pre = new point;
		pre->vec = new short[m];
		// pre = now - A[col];
		for (i = 0; i < m; i++) { pre->vec[i] = now->vec[i] - A[col].vec[i]; }
		pre->uniqueID = now->uniqueID - A[col].uniqueID;

		feasible = true;
		for (i = 0; i < m; i++) { if (A[col].vec[i] > now->vec[i]) { feasible = false; break; } }

		if (feasible == false)
		{
			StorePoint(now, cnt);
			// lsm[cnt].insert(now);
			nowvalue = now->valfun;
		}
		else
		{
			pre->valfun = value(pre, cnt); // z_{k}(\beta - a_k)

			nowvalue = max(pre->valfun + c[col], (int)(*it)->valfun);

			now->valfun = nowvalue;

			if ((*it)->valfun > pre->valfun + c[col] || FindPoint(pre, cnt) == true)
			{
				StorePoint(now, cnt);
				// lsm[cnt].insert(now);
			}
			else
			{
				it++;
				continue;
			}
		}

		//find new ones

		for (i = 0; i < u_k; i++)
		{
			// larger = now + A[col], nowvalue = z_{k}(\beta + (t-1)a_k)
			point* larger = new point;
			larger->vec = new short[m];
			for (j = 0; j < m; j++) { larger->vec[j] = now->vec[j] + (i + 1) * A[col].vec[j]; }
			larger->uniqueID = (now->uniqueID) + (i + 1) * (A[col].uniqueID);

			feasible = true;
			for (j = 0; j < m; j++) { if (larger->vec[j] > ub.vec[j]) { feasible = false; break; } }

			if (feasible == false)
			{
				break;
			}

			prevalue = value(larger, 1 - cnt); // z_{k-1}(\beta + t a_k)

			larger->valfun = max(nowvalue + c[col], (int)prevalue); // z_k(\beta + t a_k)

			if (prevalue < nowvalue + c[col])
			{
				StorePoint(larger, cnt);
				// lsm[cnt].insert(larger);
				nowvalue = larger->valfun;
			}
			else if (FindPoint(larger, 1 - cnt) == true)
			{
				StorePoint(larger, cnt);
				// lsm[cnt].insert(larger);
				nowvalue = larger->valfun;
			}
			else
			{
				break;
			}
		}

		it++;
	}
}

void readdata(string datafile_name)
{
	int i, j;
	short* cc = NULL;

	ifstream datafile;
	datafile.open(datafile_name, ifstream::in);
	if (!datafile.is_open()) {
		cout << "Problem: data file is not available!" << endl;
		exit(1);
	}

	// m: number of rows
	// n: number of columns

	datafile >> m >> n;

	A.resize(n);

	maxrhs = new short[m];
	for (i = 0; i < m; i++) { datafile >> maxrhs[i]; }

	ub.vec = new short[m];
	for (i = 0; i < m; i++) { ub.vec[i] = maxrhs[i]; }

	indinc = new double[m];
	indinc[m - 1] = 1;
	for (i = m - 2; i >= 0; i--) { indinc[i] = indinc[i + 1] * (maxrhs[i + 1] + 1); }

	cc = new short[n];
	for (j = 0; j < n; j++) { datafile >> cc[j]; }

	vector<vector<int>> AA(n, vector<int>(m, 0));

	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) { datafile >> AA[j][i]; }
	}
	datafile.close();

	/*for (j = 0; j < n; j++) cout << cc[j] << " ";
	cout << endl;

	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++)
			cout << AA[j][i] << " ";
		cout << endl;
	}
	cout << endl;*/

	/////////////////////////////////////////////////////////////////////
	// Keep the column order unchanged
	/*
	for (j = 0; j < n; j++) {
		A[j].vec = new short[m];
		for (i = 0; i < m; i++)
			A[j].vec[i] = AA[j][i];
	}

	for (j = 0; j < n; j++) { A[j].uniqueID = src(&A[j]); }

	c = new short[n];
	for (j = 0; j < n; j++) { c[j] = cc[j]; }
	*/
	/////////////////////////////////////////////////////////////////////


	/////////////////////////////////////////////////////////////////////
	// Sort the columns in decreasing order of u_k
	/*
	vector<int> u_k(n, INFTY); // the maximum value that each variable can take with the max rhs

	for (j = 0; j < n; j++)
	{
		for (i = 0; i < m; i++)
		{
			if (AA[j][i] == 0) { continue; }

			int tempint = floor(maxrhs[i] / AA[j][i]);

			if (tempint < u_k[j]) { u_k[j] = tempint; }
		}
	}

	// Vector to store element with respective present index
	vector<pair<int, int>> vp;

	// Inserting element in pair vector to keep track of previous indexes
	for (j = 0; j < n; j++) {
		vp.push_back(make_pair(u_k[j], j));
	}

	// Sorting pair vector
	sort(vp.begin(), vp.end());

	for (j = 0; j < n; j++) {
		A[j].vec = new short[m];
		for (i = 0; i < m; i++)
			A[j].vec[i] = AA[vp[n - 1 - j].second][i];
	}

	for (j = 0; j < n; j++) { A[j].uniqueID = src(&A[j]); }

	c = new short[n];
	for (j = 0; j < n; j++) { c[j] = cc[vp[n - 1 - j].second]; }
	*/
	/////////////////////////////////////////////////////////////////////


	/////////////////////////////////////////////////////////////////////
	// Sort the columns in increasing order of u_k
	/*
	vector<int> u_k(n, INFTY); // the maximum value that each variable can take with the max rhs

	for (j = 0; j < n; j++) 
	{
		for (i = 0; i < m; i++)
		{
			if (AA[j][i] == 0) { continue; }

			int tempint = floor(maxrhs[i] / AA[j][i]);

			if (tempint < u_k[j]) { u_k[j] = tempint; }
		}
	}
	
	// Vector to store element with respective present index
	vector<pair<int, int>> vp;

	// Inserting element in pair vector to keep track of previous indexes
	for (j = 0; j < n; j++) {
		vp.push_back(make_pair(u_k[j], j));
	}

	// Sorting pair vector
	sort(vp.begin(), vp.end());

	for (j = 0; j < n; j++) {
		A[j].vec = new short[m];
		for (i = 0; i < m; i++)
			A[j].vec[i] = AA[vp[j].second][i];
	}

	for (j = 0; j < n; j++) { A[j].uniqueID = src(&A[j]); }

	c = new short[n];
	for (j = 0; j < n; j++) { c[j] = cc[vp[j].second]; }
	*/
	/////////////////////////////////////////////////////////////////////


	/////////////////////////////////////////////////////////////////////
	// Sort the columns in lexicographical order
	
	// Vector to store element with respective present index
	vector<pair<vector<int>, int>> vp;

	// Inserting element in pair vector to keep track of previous indexes
	for (j = 0; j < n; j++) {
		vp.push_back(make_pair(AA[j], j));
	}

	// Sorting pair vector
	sort(vp.begin(), vp.end());

	for (j = 0; j < n; j++) {
		A[j].vec = new short[m];
		for (i = 0; i < m; i++)
			A[j].vec[i] = vp[j].first[i];
	}

	for (j = 0; j < n; j++) { A[j].uniqueID = src(&A[j]); }

	c = new short[n];
	for (j = 0; j < n; j++) { c[j] = cc[vp[j].second]; }
	
	/////////////////////////////////////////////////////////////////////


	/////////////////////////////////////////////////////////////////////
	// Sort the columns in reverse lexicographical order
	/*
	// Vector to store element with respective present index
	vector<pair<vector<int>, int>> vp;

	// Inserting element in pair vector to keep track of previous indexes
	for (j = 0; j < n; j++) {
		vp.push_back(make_pair(AA[j], j));
	}

	// Sorting pair vector
	sort(vp.begin(), vp.end());

	for (j = 0; j < n; j++) {
		A[j].vec = new short[m];
		for (i = 0; i < m; i++)
			A[j].vec[i] = vp[n-1-j].first[i];
	}

	for (j = 0; j < n; j++) { A[j].uniqueID = src(&A[j]); }

	c = new short[n];
	for (j = 0; j < n; j++) { c[j] = cc[vp[n - 1 - j].second]; }
	*/
	/////////////////////////////////////////////////////////////////////

	/*for (j = 0; j < n; j++) cout << c[j] << " ";
	cout << endl << endl;

	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++)
			cout << A[j].vec[i] << " ";
		cout << endl;
	}
	cout << endl;*/

	delete[] cc; cc = NULL;
}


int main(int argc, char* argv[])
{
	int i, k;
	int status = 0;
	int t;
	int deletedcolumns;
	bool feasible;

	int* indices;
	double* newrhs;

	vector<point*>::iterator it;

	time_t optstart, optend;

	//////////////////////////////////////////////////////////////////////////
	string folder = "E:\\Data\\Andy instances\\Instances_copy\\";
	// string instance = "IC-T4-3-second-stage.txt";
	string instance = argv[1];
	string instancedir = folder + instance;
	/*string resultdir = folder + "Result_cbc_unordered.txt";
	string exceldir = folder + "Result_cbc_unordered_excel.txt";*/
	string resultdir = folder + "Result_cbc_lexicographical.txt";
	string exceldir = folder + "Result_cbc_lexicographical_excel.txt";
	/*string resultdir = folder + "Result_cbc_reverse_lexicographical.txt";
	string exceldir = folder + "Result_cbc_reverse_lexicographical_excel.txt";*/
	/*string resultdir = folder + "Result_cbc_increasing_u_k.txt";
	string exceldir = folder + "Result_cbc_increasing_u_k_excel.txt";*/
	/*string resultdir = folder + "Result_cbc_decreasing_u_k.txt";
	string exceldir = folder + "Result_cbc_decreasing_u_k_excel.txt";*/
	//////////////////////////////////////////////////////////////////////////

	ofstream resultfile, excelfile;
	resultfile.open(resultdir, ofstream::out | ofstream::app);
	excelfile.open(exceldir, ofstream::out | ofstream::app);

	tlimit2 = atof(argv[2]);

	resultfile << "Problem: " << instance << endl << endl;
	cout << "Problem: " << instance << endl << endl;

	readdata(instancedir); //// The order of columns are set in this function ////

	indices = new int[n];
	newrhs = new double[m];

	//////////////////////////////////////////////////////////////////////////
	// column by column begins
	// the first column

	optstart = clock(); // start recording running time

	cout << "  The number of columns: " << n << endl;
	resultfile << "  The number of columns: " << n << endl;

	deletedcolumns = 0;

	t = 0;
	
	while (1)
	{
		point* p = new point;
		p->vec = new short[m];

		for (i = 0; i < m; i++) { p->vec[i] = t * A[0].vec[i]; }

		feasible = true;
		for (i = 0; i < m; i++) { if (p->vec[i] > ub.vec[i]) { feasible = false; break; } }

		if (feasible == false) { break; }

		p->valfun = t * c[0];
		p->uniqueID = t * A[0].uniqueID;
		lsm[cnt].push_back(p);
		t++;
	}

	// the rest columns

	for (k = 1; k < n; k++)
	{
		feasible = true;
		for (i = 0; i < m; i++) { if (A[k].vec[i] > ub.vec[i]) { feasible = false; break; } }

		if (feasible == false || value(&A[k], cnt) >= c[k])
		{ 
			deletedcolumns++;
			continue; 
		}

		cnt = 1 - cnt; 

		lsm[cnt].clear();

		update(k);
	}

	optend = clock();

	cout << "  The number of columns deleted: " << deletedcolumns << endl;
	cout << "  The number of remaining columns: " << n - deletedcolumns << endl << endl;

	cout << "  Number of lsm: " << lsm[cnt].size() << endl;
	cout << "  Total running time: " << (double)(optend - optstart) / CLOCKS_PER_SEC << endl << endl;

	resultfile << "  The number of columns deleted: " << deletedcolumns << endl;
	resultfile << "  The number of remaining columns: " << n - deletedcolumns << endl << endl;

	resultfile << "  Number of lsm: " << lsm[cnt].size() << endl;
	resultfile << "  Total running time: " << (double)(optend - optstart) / CLOCKS_PER_SEC << endl << endl;

	excelfile << deletedcolumns << ' ';
	excelfile << lsm[cnt].size() << ' ';
	excelfile << (double)(optend - optstart) / CLOCKS_PER_SEC << ' ';

	excelfile << endl;

	//////////////////////////////////////////////////////////////////////////


	resultfile.close();
	excelfile.close();

	lsm[0].clear();
	lsm[1].clear();

	delete[] c; c = NULL;
	delete[] maxrhs; maxrhs = NULL;
	delete[] indinc; indinc = NULL;
	delete[] indices; indices = NULL;
	delete[] newrhs; newrhs = NULL;

	return 0;
}