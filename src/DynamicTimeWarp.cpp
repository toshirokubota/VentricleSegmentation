#include <DynamicTimeWarp.h>

float
DynamicTimeWarp(vector<float>& x, vector<float>& y, vector<float>& params)
{
	vector<vector<float>> C(x.size() + 1); //cost table
	vector<vector<int>> J(x.size() + 1); //predecessor table
	for (int i = 0; i < C.size(); ++i)
	{
		C[i] = vector<float>(y.size() + 1, std::numeric_limits<float>::infinity());
		J[i] = vector<int>(y.size() + 1, -1);
	}
	C[0][0] = 0;
	for (int i = 0; i < x.size(); ++i)
	{
		for (int j = 0; j < y.size(); ++j)
		{
			float df = x[i] - y[j];
			float cost = abs(df);
			if (C[i][j] <= C[i + 1][j] && C[i][j] <= C[i][j + 1])
			{
				C[i + 1][j + 1] = C[i][j] + cost;
				J[i + 1][j + 1] = j;
			}
			else if (C[i + 1][j] <= C[i][j + 1]){
				C[i + 1][j + 1] = C[i + 1][j] + cost;
				J[i + 1][j + 1] = J[i + 1][j];
			}
			else {
				C[i + 1][j + 1] = C[i][j + 1] + cost;
				J[i + 1][j + 1] = J[i][j + 1];
			}
		}
	}
	return C[x.size()][y.size()];
}
