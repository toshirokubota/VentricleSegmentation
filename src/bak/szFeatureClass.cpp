#include <szFeatureClass.h>

ostream& operator << (ostream& out, const CFeature& f)
{
	out << f.Name << "\t";
	out << f.Value << endl;

	return out;
}

#ifdef MEX_DLL
mxArray*
StoreFeatures(const vector<CFeature>& features) {
	mxArray* plhs;
	const int dims[2]={1, features.size()};
	const int NumberOfFields = 2;
	const char* fieldNames[] = {"name", "value"};
	/* Create a 1-by-n array of structs. */ 
	plhs = mxCreateStructArray(2, dims, NumberOfFields, fieldNames);
	int name_field = mxGetFieldNumber(plhs,"name");
	int value_field = mxGetFieldNumber(plhs,"value");
	for (int i=0; i<features.size(); i++) {
		mxArray *field_value;
		mxSetFieldByNumber(plhs, i, name_field, mxCreateString(features[i].Name.c_str()));
		field_value = mxCreateDoubleMatrix(1,1,mxREAL);
		*mxGetPr(field_value) = features[i].Value;
		mxSetFieldByNumber(plhs, i, value_field, field_value);
	}

	return plhs;
}
#endif