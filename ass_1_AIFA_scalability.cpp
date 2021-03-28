// A C++ Program to implement A* Search Algorithm
#include <bits/stdc++.h>
#include <cstdlib>
#include <chrono>
using namespace std;
using namespace std::chrono;

/************
 *************
 CHANGE THESE PARAMS  FOR SCALABILITY TESTING
 **************
 *************
 */

//#define ROW 100
//#define COL 200
//#define N 10
//#define MAX 1000 

#define ROW 200
#define COL 400
#define N 20

#define verbose 0
#define INF 0x7FFFFFFF

// Creating a shortcut for int, int pair type
typedef pair<int, int> Pair;

// Creating a shortcut for pair<int, pair<int, int>> type
typedef pair<double, pair<int, int> > pPair;

// A structure to hold the neccesary parameters
struct cell {
	// Row and Column index of its parent
	// Note that 0 <= i <= ROW-1 & 0 <= j <= COL-1
	int parent_i, parent_j;
	// f = g + h
	double f, g, h;
};

// A Utility Function to check whether given cell (row, col)
// is a valid cell or not.
bool isValid(int row, int col)
{
	// Returns true if row number and column number
	// is in range
	return (row >= 0) && (row < ROW) && (col >= 0)
		&& (col < COL);
}

// A Utility Function to check whether the given cell is
// blocked or not
bool isUnBlocked(int grid[][COL], int row, int col)
{
	// Returns true if the cell is not blocked else false
	if (grid[row][col] == 1)
		return (true);
	else
		return (false);
}

// A Utility Function to check whether destination cell has
// been reached or not
bool isDestination(int row, int col, Pair dest)
{
	if (row == dest.first && col == dest.second)
		return (true);
	else
		return (false);
}

// A Utility Function to calculate the 'h' heuristics.
double calculateHValue(int row, int col, Pair dest)
{
	// Return using the distance formula
	return ((double)sqrt(
		(row - dest.first) * (row - dest.first)
		+ (col - dest.second) * (col - dest.second)));
}

// A Utility Function to trace the path from the source
// to destination
stack<Pair> tracePath(cell cellDetails[][COL], Pair dest)
{
        //printf("\nThe Path is ");
	int row = dest.first;
	int col = dest.second;

	stack<Pair> Path;

	while (!(cellDetails[row][col].parent_i == row
			&& cellDetails[row][col].parent_j == col)) {
		Path.push(make_pair(row, col));
		int temp_row = cellDetails[row][col].parent_i;
		int temp_col = cellDetails[row][col].parent_j;
		row = temp_row;
		col = temp_col;
	}

	Path.push(make_pair(row, col));
	return Path;
	//while (!Path.empty()) {
	//	pair<int, int> p = Path.top();
	//	Path.pop();
	//	printf("-> (%d,%d) ", p.first, p.second);
	//}

}

// A Function to find the shortest path between
// a given source cell to a destination cell according
// to A* Search Algorithm
stack<Pair> aStarSearchSwa(int grid[][COL], Pair src, Pair dest)
{       
	stack<Pair> empty;
	// If the source is out of range
	if (isValid(src.first, src.second) == false) {
		printf("Source is invalid\n");
		return empty;
	}

	// If the destination is out of range
	if (isValid(dest.first, dest.second) == false) {
		printf("Destination is invalid\n");
		return empty;
	}

	// Either the source or the destination is blocked
	if (isUnBlocked(grid, src.first, src.second) == false
		|| isUnBlocked(grid, dest.first, dest.second)
			== false) {
		printf("Source or the destination is blocked\n");
		return empty;
	}

	// If the destination cell is the same as source cell
	if (isDestination(src.first, src.second, dest)
		== true) {
		printf("We are already at the destination\n");
		return empty;
	}

	// Create a closed list and initialise it to false which
	// means that no cell has been included yet This closed
	// list is implemented as a boolean 2D array
	bool closedList[ROW][COL];
	memset(closedList, false, sizeof(closedList));

	// Declare a 2D array of structure to hold the details
	// of that cell
	cell cellDetails[ROW][COL];

	int i, j;

	for (i = 0; i < ROW; i++) {
		for (j = 0; j < COL; j++) {
			cellDetails[i][j].f = FLT_MAX;
			cellDetails[i][j].g = FLT_MAX;
			cellDetails[i][j].h = FLT_MAX;
			cellDetails[i][j].parent_i = -1;
			cellDetails[i][j].parent_j = -1;
		}
	}

	// Initialising the parameters of the starting node
	i = src.first, j = src.second;
	cellDetails[i][j].f = 0.0;
	cellDetails[i][j].g = 0.0;
	cellDetails[i][j].h = 0.0;
	cellDetails[i][j].parent_i = i;
	cellDetails[i][j].parent_j = j;

	/*
	Create an open list having information as-
	<f, <i, j>>
	where f = g + h,
	and i, j are the row and column index of that cell
	Note that 0 <= i <= ROW-1 & 0 <= j <= COL-1
	This open list is implenented as a set of pair of
	pair.*/
	set<pPair> openList;

	// Put the starting cell on the open list and set its
	// 'f' as 0
	openList.insert(make_pair(0.0, make_pair(i, j)));

	// We set this boolean value as false as initially
	// the destination is not reached.
	bool foundDest = false;

	while (!openList.empty()) {
		pPair p = *openList.begin();

		// Remove this vertex from the open list
		openList.erase(openList.begin());

		// Add this vertex to the closed list
		i = p.second.first;
		j = p.second.second;
		closedList[i][j] = true;

		/*
		Generating all the 8 successor of this cell

			N.W N N.E
			\ | /
				\ | /
			W----Cell----E
				/ | \
				/ | \
			S.W S S.E

		Cell-->Popped Cell (i, j)
		N --> North	 (i-1, j)
		S --> South	 (i+1, j)
		E --> East	 (i, j+1)
		W --> West		 (i, j-1)
		N.E--> North-East (i-1, j+1)
		N.W--> North-West (i-1, j-1)
		S.E--> South-East (i+1, j+1)
		S.W--> South-West (i+1, j-1)*/

		// To store the 'g', 'h' and 'f' of the 8 successors
		double gNew, hNew, fNew;

		//----------- 1st Successor (North) ------------

		// Only process this cell if this is a valid one
		if (isValid(i - 1, j) == true) {
			// If the destination cell is the same as the
			// current successor
			if (isDestination(i - 1, j, dest) == true) {
				// Set the Parent of the destination cell
				cellDetails[i - 1][j].parent_i = i;
				cellDetails[i - 1][j].parent_j = j;
				//printf("The destination cell is found\n");
				foundDest = true;
				return tracePath(cellDetails, dest);
			}
			// If the successor is already on the closed
			// list or if it is blocked, then ignore it.
			// Else do the following
			else if (closedList[i - 1][j] == false
					&& isUnBlocked(grid, i - 1, j)
							== true) {
				gNew = cellDetails[i][j].g + 1.0;
				hNew = calculateHValue(i - 1, j, dest);
				fNew = gNew + hNew;

				// If it isn’t on the open list, add it to
				// the open list. Make the current square
				// the parent of this square. Record the
				// f, g, and h costs of the square cell
				//			 OR
				// If it is on the open list already, check
				// to see if this path to that square is
				// better, using 'f' cost as the measure.
				if (cellDetails[i - 1][j].f == FLT_MAX
					|| cellDetails[i - 1][j].f > fNew) {
					openList.insert(make_pair(
						fNew, make_pair(i - 1, j)));

					// Update the details of this cell
					cellDetails[i - 1][j].f = fNew;
					cellDetails[i - 1][j].g = gNew;
					cellDetails[i - 1][j].h = hNew;
					cellDetails[i - 1][j].parent_i = i;
					cellDetails[i - 1][j].parent_j = j;
				}
			}
		}

		//----------- 2nd Successor (South) ------------

		// Only process this cell if this is a valid one
		if (isValid(i + 1, j) == true) {
			// If the destination cell is the same as the
			// current successor
			if (isDestination(i + 1, j, dest) == true) {
				// Set the Parent of the destination cell
				cellDetails[i + 1][j].parent_i = i;
				cellDetails[i + 1][j].parent_j = j;
				//printf("The destination cell is found\n");
				foundDest = true;
                                return tracePath(cellDetails, dest);
			}
			// If the successor is already on the closed
			// list or if it is blocked, then ignore it.
			// Else do the following
			else if (closedList[i + 1][j] == false
					&& isUnBlocked(grid, i + 1, j)
							== true) {
				gNew = cellDetails[i][j].g + 1.0;
				hNew = calculateHValue(i + 1, j, dest);
				fNew = gNew + hNew;

				// If it isn’t on the open list, add it to
				// the open list. Make the current square
				// the parent of this square. Record the
				// f, g, and h costs of the square cell
				//			 OR
				// If it is on the open list already, check
				// to see if this path to that square is
				// better, using 'f' cost as the measure.
				if (cellDetails[i + 1][j].f == FLT_MAX
					|| cellDetails[i + 1][j].f > fNew) {
					openList.insert(make_pair(
						fNew, make_pair(i + 1, j)));
					// Update the details of this cell
					cellDetails[i + 1][j].f = fNew;
					cellDetails[i + 1][j].g = gNew;
					cellDetails[i + 1][j].h = hNew;
					cellDetails[i + 1][j].parent_i = i;
					cellDetails[i + 1][j].parent_j = j;
				}
			}
		}

		//----------- 3rd Successor (East) ------------

		// Only process this cell if this is a valid one
		if (isValid(i, j + 1) == true) {
			// If the destination cell is the same as the
			// current successor
			if (isDestination(i, j + 1, dest) == true) {
				// Set the Parent of the destination cell
				cellDetails[i][j + 1].parent_i = i;
				cellDetails[i][j + 1].parent_j = j;
				//printf("The destination cell is found\n");
				foundDest = true;
                                return tracePath(cellDetails, dest);
			}

			// If the successor is already on the closed
			// list or if it is blocked, then ignore it.
			// Else do the following
			else if (closedList[i][j + 1] == false
					&& isUnBlocked(grid, i, j + 1)
							== true) {
				gNew = cellDetails[i][j].g + 1.0;
				hNew = calculateHValue(i, j + 1, dest);
				fNew = gNew + hNew;

				// If it isn’t on the open list, add it to
				// the open list. Make the current square
				// the parent of this square. Record the
				// f, g, and h costs of the square cell
				//			 OR
				// If it is on the open list already, check
				// to see if this path to that square is
				// better, using 'f' cost as the measure.
				if (cellDetails[i][j + 1].f == FLT_MAX
					|| cellDetails[i][j + 1].f > fNew) {
					openList.insert(make_pair(
						fNew, make_pair(i, j + 1)));

					// Update the details of this cell
					cellDetails[i][j + 1].f = fNew;
					cellDetails[i][j + 1].g = gNew;
					cellDetails[i][j + 1].h = hNew;
					cellDetails[i][j + 1].parent_i = i;
					cellDetails[i][j + 1].parent_j = j;
				}
			}
		}

		//----------- 4th Successor (West) ------------

		// Only process this cell if this is a valid one
		if (isValid(i, j - 1) == true) {
			// If the destination cell is the same as the
			// current successor
			if (isDestination(i, j - 1, dest) == true) {
				// Set the Parent of the destination cell
				cellDetails[i][j - 1].parent_i = i;
				cellDetails[i][j - 1].parent_j = j;
				//printf("The destination cell is found\n");
				foundDest = true;
                                return tracePath(cellDetails, dest);
			}

			// If the successor is already on the closed
			// list or if it is blocked, then ignore it.
			// Else do the following
			else if (closedList[i][j - 1] == false
					&& isUnBlocked(grid, i, j - 1)
							== true) {
				gNew = cellDetails[i][j].g + 1.0;
				hNew = calculateHValue(i, j - 1, dest);
				fNew = gNew + hNew;

				// If it isn’t on the open list, add it to
				// the open list. Make the current square
				// the parent of this square. Record the
				// f, g, and h costs of the square cell
				//			 OR
				// If it is on the open list already, check
				// to see if this path to that square is
				// better, using 'f' cost as the measure.
				if (cellDetails[i][j - 1].f == FLT_MAX
					|| cellDetails[i][j - 1].f > fNew) {
					openList.insert(make_pair(
						fNew, make_pair(i, j - 1)));

					// Update the details of this cell
					cellDetails[i][j - 1].f = fNew;
					cellDetails[i][j - 1].g = gNew;
					cellDetails[i][j - 1].h = hNew;
					cellDetails[i][j - 1].parent_i = i;
					cellDetails[i][j - 1].parent_j = j;
				}
			}
		}

		//----------- 5th Successor (North-East)
		//------------

		// Only process this cell if this is a valid one
		if (isValid(i - 1, j + 1) == true) {
			// If the destination cell is the same as the
			// current successor
			if (isDestination(i - 1, j + 1, dest) == true) {
				// Set the Parent of the destination cell
				cellDetails[i - 1][j + 1].parent_i = i;
				cellDetails[i - 1][j + 1].parent_j = j;
				//printf("The destination cell is found\n");
				foundDest = true;
                                return tracePath(cellDetails, dest);
			}

			// If the successor is already on the closed
			// list or if it is blocked, then ignore it.
			// Else do the following
			else if (closedList[i - 1][j + 1] == false
					&& isUnBlocked(grid, i - 1, j + 1)
							== true) {
				gNew = cellDetails[i][j].g + 1.414;
				hNew = calculateHValue(i - 1, j + 1, dest);
				fNew = gNew + hNew;

				// If it isn’t on the open list, add it to
				// the open list. Make the current square
				// the parent of this square. Record the
				// f, g, and h costs of the square cell
				//			 OR
				// If it is on the open list already, check
				// to see if this path to that square is
				// better, using 'f' cost as the measure.
				if (cellDetails[i - 1][j + 1].f == FLT_MAX
					|| cellDetails[i - 1][j + 1].f > fNew) {
					openList.insert(make_pair(
						fNew, make_pair(i - 1, j + 1)));

					// Update the details of this cell
					cellDetails[i - 1][j + 1].f = fNew;
					cellDetails[i - 1][j + 1].g = gNew;
					cellDetails[i - 1][j + 1].h = hNew;
					cellDetails[i - 1][j + 1].parent_i = i;
					cellDetails[i - 1][j + 1].parent_j = j;
				}
			}
		}

		//----------- 6th Successor (North-West)
		//------------

		// Only process this cell if this is a valid one
		if (isValid(i - 1, j - 1) == true) {
			// If the destination cell is the same as the
			// current successor
			if (isDestination(i - 1, j - 1, dest) == true) {
				// Set the Parent of the destination cell
				cellDetails[i - 1][j - 1].parent_i = i;
				cellDetails[i - 1][j - 1].parent_j = j;
				//printf("The destination cell is found\n");
				foundDest = true;
                                return tracePath(cellDetails, dest);
			}

			// If the successor is already on the closed
			// list or if it is blocked, then ignore it.
			// Else do the following
			else if (closedList[i - 1][j - 1] == false
					&& isUnBlocked(grid, i - 1, j - 1)
							== true) {
				gNew = cellDetails[i][j].g + 1.414;
				hNew = calculateHValue(i - 1, j - 1, dest);
				fNew = gNew + hNew;

				// If it isn’t on the open list, add it to
				// the open list. Make the current square
				// the parent of this square. Record the
				// f, g, and h costs of the square cell
				//			 OR
				// If it is on the open list already, check
				// to see if this path to that square is
				// better, using 'f' cost as the measure.
				if (cellDetails[i - 1][j - 1].f == FLT_MAX
					|| cellDetails[i - 1][j - 1].f > fNew) {
					openList.insert(make_pair(
						fNew, make_pair(i - 1, j - 1)));
					// Update the details of this cell
					cellDetails[i - 1][j - 1].f = fNew;
					cellDetails[i - 1][j - 1].g = gNew;
					cellDetails[i - 1][j - 1].h = hNew;
					cellDetails[i - 1][j - 1].parent_i = i;
					cellDetails[i - 1][j - 1].parent_j = j;
				}
			}
		}

		//----------- 7th Successor (South-East)
		//------------

		// Only process this cell if this is a valid one
		if (isValid(i + 1, j + 1) == true) {
			// If the destination cell is the same as the
			// current successor
			if (isDestination(i + 1, j + 1, dest) == true) {
				// Set the Parent of the destination cell
				cellDetails[i + 1][j + 1].parent_i = i;
				cellDetails[i + 1][j + 1].parent_j = j;
				//printf("The destination cell is found\n");
				foundDest = true;
                                return tracePath(cellDetails, dest);
			}

			// If the successor is already on the closed
			// list or if it is blocked, then ignore it.
			// Else do the following
			else if (closedList[i + 1][j + 1] == false
					&& isUnBlocked(grid, i + 1, j + 1)
							== true) {
				gNew = cellDetails[i][j].g + 1.414;
				hNew = calculateHValue(i + 1, j + 1, dest);
				fNew = gNew + hNew;

				// If it isn’t on the open list, add it to
				// the open list. Make the current square
				// the parent of this square. Record the
				// f, g, and h costs of the square cell
				//			 OR
				// If it is on the open list already, check
				// to see if this path to that square is
				// better, using 'f' cost as the measure.
				if (cellDetails[i + 1][j + 1].f == FLT_MAX
					|| cellDetails[i + 1][j + 1].f > fNew) {
					openList.insert(make_pair(
						fNew, make_pair(i + 1, j + 1)));

					// Update the details of this cell
					cellDetails[i + 1][j + 1].f = fNew;
					cellDetails[i + 1][j + 1].g = gNew;
					cellDetails[i + 1][j + 1].h = hNew;
					cellDetails[i + 1][j + 1].parent_i = i;
					cellDetails[i + 1][j + 1].parent_j = j;
				}
			}
		}

		//----------- 8th Successor (South-West)
		//------------

		// Only process this cell if this is a valid one
		if (isValid(i + 1, j - 1) == true) {
			// If the destination cell is the same as the
			// current successor
			if (isDestination(i + 1, j - 1, dest) == true) {
				// Set the Parent of the destination cell
				cellDetails[i + 1][j - 1].parent_i = i;
				cellDetails[i + 1][j - 1].parent_j = j;
				//printf("The destination cell is found\n");
				foundDest = true;
                                return tracePath(cellDetails, dest);
			}

			// If the successor is already on the closed
			// list or if it is blocked, then ignore it.
			// Else do the following
			else if (closedList[i + 1][j - 1] == false
					&& isUnBlocked(grid, i + 1, j - 1)
							== true) {
				gNew = cellDetails[i][j].g + 1.414;
				hNew = calculateHValue(i + 1, j - 1, dest);
				fNew = gNew + hNew;

				// If it isn’t on the open list, add it to
				// the open list. Make the current square
				// the parent of this square. Record the
				// f, g, and h costs of the square cell
				//			 OR
				// If it is on the open list already, check
				// to see if this path to that square is
				// better, using 'f' cost as the measure.
				if (cellDetails[i + 1][j - 1].f == FLT_MAX
					|| cellDetails[i + 1][j - 1].f > fNew) {
					openList.insert(make_pair(
						fNew, make_pair(i + 1, j - 1)));

					// Update the details of this cell
					cellDetails[i + 1][j - 1].f = fNew;
					cellDetails[i + 1][j - 1].g = gNew;
					cellDetails[i + 1][j - 1].h = hNew;
					cellDetails[i + 1][j - 1].parent_i = i;
					cellDetails[i + 1][j - 1].parent_j = j;
				}
			}
		}
	}

	// When the destination cell is not found and the open
	// list is empty, then we conclude that we failed to
	// reach the destiantion cell. This may happen when the
	// there is no way to destination cell (due to
	// blockages)
	if (foundDest == false)
		printf("Failed to find the Destination Cell\n");

	return empty;
}


// state space tree node
struct Node
{
	// stores parent node of current node
	// helps in tracing path when answer is found
	Node* parent;

	// contains cost for ancestors nodes
	// including current node
	int pathCost;

	// contains least promising cost
	int cost;

	// contain worker number
	int workerID;

	// contains Job ID
	int jobID;

	// Boolean array assigned will contains
	// info about available jobs
	bool assigned[N];
};

// Function to allocate a new search tree node
// Here Person x is assigned to job y
Node* newNode(int x, int y, bool assigned[],
			Node* parent)
{
	Node* node = new Node;

	for (int j = 0; j < N; j++)
		node->assigned[j] = assigned[j];
	node->assigned[y] = true;

	node->parent = parent;
	node->workerID = x;
	node->jobID = y;

	return node;
}

// Function to calculate the least promising cost
// of node after worker x is assigned to job y.
int calculateCost(int costMatrix[N][N], int x,
				int y, bool assigned[])
{
	int cost = 0;

	// to store unavailable jobs
	bool available[N] = {true};

	// start from next worker
	for (int i = x + 1; i < N; i++)
	{
		int min = INT_MAX, minIndex = -1;

		// do for each job
		for (int j = 0; j < N; j++)
		{
			// if job is unassigned
			if (!assigned[j] && available[j] &&
				costMatrix[i][j] < min)
			{
				// store job number
				minIndex = j;

				// store cost
				min = costMatrix[i][j];
			}
		}

		// add cost of next worker
		cost += min;

		// job becomes unavailable
		available[minIndex] = false;
	}

	return cost;
}

// Comparison object to be used to order the heap
struct comp
{
	bool operator()(const Node* lhs,
				const Node* rhs) const
	{
		return lhs->cost > rhs->cost;
	}
};

// print Assignments
void printAssignments(Node *min)
{
	if(min->parent==NULL)
		return;

	printAssignments(min->parent);
	//cout << "Assign Robot " << char(min->workerID + 'A')
	cout << endl;
	cout << "Assign Robot-" << min->workerID 
		<< " to Task-" << min->jobID << endl;

}

// Finds minimum cost using Branch and Bound.
int findMinCost(int costMatrix[N][N])
{
	// Create a priority queue to store live nodes of
	// search tree;
	priority_queue<Node*, std::vector<Node*>, comp> pq;

	// initailize heap to dummy node with cost 0
	bool assigned[N] = {false};
	Node* root = newNode(-1, -1, assigned, NULL);
	root->pathCost = root->cost = 0;
	root->workerID = -1;

	// Add dummy node to list of live nodes;
	pq.push(root);

	// Finds a live node with least cost,
	// add its childrens to list of live nodes and
	// finally deletes it from the list.
	while (!pq.empty())
	{
	// Find a live node with least estimated cost
	Node* min = pq.top();

	// The found node is deleted from the list of
	// live nodes
	pq.pop();

	// i stores next worker
	int i = min->workerID + 1;

	// if all workers are assigned a job
	if (i == N)
	{
		printAssignments(min);
		return min->cost;
	}

	// do for each job
	for (int j = 0; j < N; j++)
	{
		// If unassigned
		if (!min->assigned[j])
		{
		// create a new tree node
		Node* child = newNode(i, j, min->assigned, min);

		// cost for ancestors nodes including current node
		child->pathCost = min->pathCost + costMatrix[i][j];

		// calculate its lower bound
		child->cost = child->pathCost +
			calculateCost(costMatrix, i, j, child->assigned);

		// Add child to list of live nodes;
		pq.push(child);
		}
	}
	}
        return -1;
}

char Result[N][N];  // used as boolean
void hungarian(int Array[N][N])
{
    int i,j;
    
    unsigned int m=N,n=N;
    int k;
    int l;
    int s;
    int col_mate[N]={0};
    int row_mate[N]={0};
    int parent_row[N]={0};
    int unchosen_row[N]={0};
    int t;
    int q;
    int row_dec[N]={0};
    int col_inc[N]={0};
    int slack[N]={0};
    int slack_row[N]={0};
    int unmatched;
    int cost=0;
    
    for (i=0;i<N;++i)
      for (j=0;j<N;++j)
        Result[i][j]=false;
    
    // Begin subtract column minima in order to start with lots of zeroes 12
    printf("Using heuristic\n");
    for (l=0;l<n;l++)
    {
      s=Array[0][l];
      for (k=1;k<n;k++)
        if (Array[k][l]<s)
          s=Array[k][l];
      cost+=s;
      if (s!=0)
        for (k=0;k<n;k++)
          Array[k][l]-=s;
    }
    // End subtract column minima in order to start with lots of zeroes 12
    
    // Begin initial state 16
    t=0;
    for (l=0;l<n;l++)
    {
      row_mate[l]= -1;
      parent_row[l]= -1;
      col_inc[l]=0;
      slack[l]=INF;
    }
    for (k=0;k<m;k++)
    {
      s=Array[k][0];
      for (l=1;l<n;l++)
        if (Array[k][l]<s)
          s=Array[k][l];
      row_dec[k]=s;
      for (l=0;l<n;l++)
        if (s==Array[k][l] && row_mate[l]<0)
        {
          col_mate[k]=l;
          row_mate[l]=k;
          if (verbose)
            printf("matching col %d==row %d\n",l,k);
          goto row_done;
        }
      col_mate[k]= -1;
      if (verbose)
        printf("node %d: unmatched row %d\n",t,k);
      unchosen_row[t++]=k;
    row_done:
      ;
    }
    // End initial state 16
     
    // Begin Hungarian algorithm 18
    if (t==0)
      goto done;
    unmatched=t;
    while (1)
    {
      if (verbose)
        printf("Matched %d rows.\n",m-t);
      q=0;
      while (1)
      {
        while (q<t)
        {
          // Begin explore node q of the forest 19
          {
            k=unchosen_row[q];
            s=row_dec[k];
            for (l=0;l<n;l++)
              if (slack[l])
              {
                int del;
                del=Array[k][l]-s+col_inc[l];
                if (del<slack[l])
                {
                  if (del==0)
                  {
                    if (row_mate[l]<0)
                      goto breakthru;
                    slack[l]=0;
                    parent_row[l]=k;
                    if (verbose)
                      printf("node %d: row %d==col %d--row %d\n",
                        t,row_mate[l],l,k);
                    unchosen_row[t++]=row_mate[l];
                  }
                  else
                  {
                    slack[l]=del;
                    slack_row[l]=k;
                  }
              }
            }
          }
          // End explore node q of the forest 19
          q++;
        }
     
        // Begin introduce a new zero into the matrix 21
        s=INF;
        for (l=0;l<n;l++)
          if (slack[l] && slack[l]<s)
            s=slack[l];
        for (q=0;q<t;q++)
          row_dec[unchosen_row[q]]+=s;
        for (l=0;l<n;l++)
          if (slack[l])
          {
            slack[l]-=s;
            if (slack[l]==0)
            {
              // Begin look at a new zero 22
              k=slack_row[l];
              if (verbose)
                printf(
                  "Decreasing uncovered elements by %d produces zero at [%d,%d]\n",
                  s,k,l);
              if (row_mate[l]<0)
              {
                for (j=l+1;j<n;j++)
                  if (slack[j]==0)
                    col_inc[j]+=s;
                goto breakthru;
              }
              else
              {
                parent_row[l]=k;
                if (verbose)
                  printf("node %d: row %d==col %d--row %d\n",t,row_mate[l],l,k);
                unchosen_row[t++]=row_mate[l];
              }
              // End look at a new zero 22
            }
          }
          else
            col_inc[l]+=s;
        // End introduce a new zero into the matrix 21
      }
    breakthru:
      // Begin update the matching 20
      if (verbose)
        printf("Breakthrough at node %d of %d!\n",q,t);
      while (1)
      {
        j=col_mate[k];
        col_mate[k]=l;
        row_mate[l]=k;
        if (verbose)
          printf("rematching col %d==row %d\n",l,k);
        if (j<0)
          break;
        k=parent_row[j];
        l=j;
      }
      // End update the matching 20
      if (--unmatched==0)
        goto done;
      // Begin get ready for another stage 17
      t=0;
      for (l=0;l<n;l++)
      {
        parent_row[l]= -1;
        slack[l]=INF;
      }
      for (k=0;k<m;k++)
        if (col_mate[k]<0)
        {
          if (verbose)
            printf("node %d: unmatched row %d\n",t,k);
          unchosen_row[t++]=k;
        }
      // End get ready for another stage 17
    }
    done:
    
    // Begin doublecheck the solution 23
    for (k=0;k<m;k++)
      for (l=0;l<n;l++)
        if (Array[k][l]<row_dec[k]-col_inc[l])
          exit(0);
    for (k=0;k<m;k++)
    {
      l=col_mate[k];
      if (l<0 || Array[k][l]!=row_dec[k]-col_inc[l])
        exit(0);
    }
    k=0;
    for (l=0;l<n;l++)
      if (col_inc[l])
        k++;
    if (k>m)
      exit(0);
    // End doublecheck the solution 23
    // End Hungarian algorithm 18
    
    for (i=0;i<m;++i)
    {
      Result[i][col_mate[i]]=true;
     /*TRACE("%d - %d\n", i, col_mate[i]);*/
    }
    for (k=0;k<m;++k)
    {
      for (l=0;l<n;++l)
      {
        /*TRACE("%d ",Array[k][l]-row_dec[k]+col_inc[l]);*/
        Array[k][l]=Array[k][l]-row_dec[k]+col_inc[l];
      }
      /*TRACE("\n");*/
    }
    for (i=0;i<m;i++)
      cost+=row_dec[i];
    for (i=0;i<n;i++)
      cost-=col_inc[i];
    printf("Cost is %d\n",cost);
}


// Driver code
int main()
{

        /* Task 2: @swarndey tests for SCALABILITY 
	 * **************************************
	 *
	 * This will require :
	 * 1. a huge grid
	 * 2. Lots of robots
	 * 3. Lots of tasks
	 * 
	 */

	/* Description of the Grid-
	1--> The cell is not blocked
	0--> The cell is blocked */

        int i = 0;
	int j = 0;
        cout << "\n starting scalability test **********************\n";
	//cout << " Enter Grid x dimension\n";
        //int ROW= 0;
        //cin >> ROW;	
	//cout << " Enter Grid y dimension\n";
        //int COL= 0;
        //cin >> COL;	
	int grid1[ROW][COL];
	cout << " Create and print the Grid\n";
        for (i = 0 ; i < ROW; i++){
		for (j=0; j< COL;j++){
			if(j*i == 100 || j*i == 200 || j*i == 300 || j*i == 350 || j*i == 525 || j*i == 512 || j*i == 699 || j*i == 888 || j*i == 4000 || j*i == 1200 || j*i == 1300 || j*i == 1425 || j*i == 1508 || j*i == 1690 || j*i == 1752)
				grid1[i][j]=0;
			else
				grid1[i][j]=1;
			cout << grid1[i][j];
		}
		cout<<"\n";
	}


        /* Task 2: Define the Robot start and goal locations */

        // These arrays are indexed by the robots
	Pair starts1[N];
	Pair goals1[N];
        
	// These arrays are indexed by the tasks
	Pair picks1[N];
	Pair drops1[N];

        /* Now Dey da will have to generate the robots and the tasks 
	 * This will be done automatically
	 */
        
	//cout << " Enter Number of Robots\n";
        //int NUM_ROBO= 0;
        //cin >> NUM_ROBO;	
	//if(NUM_ROBO>MAX){
	//	cout << " Sorry, the maximum number of robots supported is: " << MAX << " You can change in the code \n";
	//	N = MAX;
	//}	
	//else
	//	N = NUM_ROBO;	
        int rnx=0; 
        int rny=0; 
	srand(1616860213);
	for (i = 0 ; i < N; i++){
		rnx = rand() % ROW;
		rny = rand() % COL;
//              cout << " rnx, rny" << rnx << "    " << rny << "\n";
	        if(isUnBlocked(grid1, rnx, rny)){
			starts1[i] = make_pair(rnx,rny);
		}
		else{
		        i--;
			continue;
		}
	}
	for (i = 0 ; i < N; i++){
		rnx = rand() % ROW;
		rny = rand() % COL;
	        if(isUnBlocked(grid1, rnx, rny)){
			goals1[i] = make_pair(rnx,rny);
		}
		else{
		        i--;
			continue;
		}
	}
	for (i = 0 ; i < N; i++){
		rnx = rand() % ROW;
		rny = rand() % COL;
	        if(isUnBlocked(grid1, rnx, rny)){
			picks1[i] = make_pair(rnx,rny);
		}
		else{
		        i--;
			continue;
		}
	}
	for (i = 0 ; i < N; i++){
		rnx = rand() % ROW;
		rny = rand() % COL;
	        if(isUnBlocked(grid1, rnx, rny)){
			drops1[i] = make_pair(rnx,rny);
		}
		else{
		        i--;
			continue;
		}
	}

	/* Okay, now @swarnava, for a robot taking up a particular task, the
	 * total cost will be start -> pick + (pick to drop + drop to goal)
	 */
	// And the temporary cost matrix here for N tasks & N agents
	// x-cordinate represents a agent/ swabot
	// y-cordinate represents a Job / swatask
	stack<Pair> tmp;
	int tmp_cost = 0;
	int costMatrix1[N][N];
       


	// For each robot
	for(i =0 ; i < N ; i++){
		// For each task
		for(j=0; j < N; j++){
		    // Find cost for robot- i and  task j pair - start to pick
	            tmp = aStarSearchSwa(grid1, starts1[i], picks1[j]);
                    tmp_cost += tmp.size();
		    // Done @sawrdey : now next part is fixed pick to drop
		    // and drop to goal for this combo
	            tmp = aStarSearchSwa(grid1, picks1[j], drops1[j]);
                    tmp_cost += tmp.size();
	            tmp = aStarSearchSwa(grid1, drops1[j], goals1[j]);
                    tmp_cost += tmp.size();
                    // Put this in the cost matrix 
		    costMatrix1[i][j] = tmp_cost;
		}
	}



        // Print the cost matrix
	cout << "\nPrint Cost matrix\n" ;
	
	// For each robot
	for(i =0 ; i < N ; i++){
		// For each task
	        cout << "\n" ;
		for(j=0; j < N; j++){
	            cout << costMatrix1[i][j] << "\t" ;
		}
	}
        auto start = high_resolution_clock::now();
	cout << "\n Optimal Search DFBB Starting ..........\n";
	cout << "\nOptimal Cost is :"  << findMinCost(costMatrix1);
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
 

        
	
	
	auto start1 = high_resolution_clock::now();
	cout << "\n Heuristic Search Starting ..........\n";
	hungarian(costMatrix1);
	auto stop1 = high_resolution_clock::now();
	auto duration1 = duration_cast<microseconds>(stop1 - start1);

        for (i=0;i<N;++i)
		for (j=0;j<N;++j){
			if (Result[i][j])
                            printf("Assign Robot-%d to  Task-%d \n",i,j);
                }  


	cout << "\n";
	cout << "\n";
	cout << "\n";
	cout << "Time taken by baranch and bound optimal search: " << duration.count() << endl;
	cout << "\n";
	cout << "\n";
	cout << "\n";

	cout << "\n";
	cout << "\n";
	cout << "\n";
	cout << "Time taken by heuristic search: " << duration1.count() << endl;
	cout << "\n";
	cout << "\n";
	cout << "\n";
	return 0;
}
