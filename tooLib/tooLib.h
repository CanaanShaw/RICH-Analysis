//
//  tooLib.h
//  RICH Analysis
//
//  Created by Canaan Shaw on 2023/6/7.
//

#ifndef tooLib_h
#define tooLib_h

#include "sharedHeader.h"

#define MAXIMUM 999

using namespace std;

// sloppy, try not to use it!!
template < typename T1, typename T2 >
vector < T1 > convert (const vector < T2 > & variable1) {

    std::stringstream ss;
    auto size = variable1.size();
    vector < T1 > result (size);

    for (int i = 0; i < size; i ++) {
        result[i] = (T1)variable1[i];
    }

    return result;
}

template < typename T1 >
vector < T1 > spinX (vector < T1 > variable1, double variable2) {

    vector < T1 > result;
    result.push_back(variable1[0]);
    result.push_back(variable1[1] * cos(variable2) - variable1[2] * sin(variable2));
    result.push_back(variable1[1] * sin(variable2) + variable1[2] * cos(variable2));
    return result;
}

template < typename T1 >
vector < T1 > spinY (vector < T1 > variable1, double variable2) {

    vector < T1 > result;
    result.push_back(variable1[0] * cos(variable2) + variable1[2] * sin(variable2));
    result.push_back(variable1[1]);
    result.push_back(-variable1[0] * sin(variable2) + variable1[2] * cos(variable2));
    return result;
}

template < typename T1 >
vector < T1 > spinZ (vector < T1 > variable1, double variable2) {

    vector < T1 > result;
    result.push_back(variable1[0] * cos(variable2) - variable1[1] * sin(variable2));
    result.push_back(variable1[0] * sin(variable2) + variable1[1] * cos(variable2));
    result.push_back(variable1[2]);
    return result;
}


template < typename T1 >
T1 square (const T1 & variable1) {

    return variable1 * variable1;
}

template < typename T1 >
T1 cube (const T1 & variable1) {

    return variable1 * variable1 * variable1;
}

template < typename T1, typename T2 >
T1 convert(const T2 & variable1) {

    std::stringstream ss;
    T1 result;

    ss << variable1;
    ss >> result;
    return result;
}

template < typename T1 >
ostream & operator << (ostream & os, const vector < T1 > & vec) {
    
    for (auto element : vec) {

        os << element << " ";
    }

    cout << endl;

    return os;
}

template < typename T1, typename T2, typename T3 >
T1 max (const T2 & variable1, const T3 & variable2) {

    double temp1 = convert <double> (variable1);
    double temp2 = convert <double> (variable2);

    if (temp1 >= temp2) {
        return convert < T1 > (temp1);
    } else {
        return convert < T1 > (temp2);
    }

    return 0;
}

template < typename T1, typename T2, typename T3 >
T1 min (const T2 & variable1, const T3 & variable2) {

    double temp1 = convert < double > (variable1);
    double temp2 = convert < double > (variable2);

    if (temp1 >= temp2) {
        return convert < T1 > (temp2);
    } else {
        return convert < T1 > (temp1);
    }

    return 0;
}

template < typename T1, typename T2 >
double distance(const vector < T1 > & variable1, const vector < T2 > & variable2) {

    double x1 = convert < double > (variable1[0]);
    double y1 = convert < double > (variable1[1]);

    double x2 = convert < double > (variable2[0]);
    double y2 = convert < double > (variable2[1]);

    return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

template < typename T1, typename T2 >
double distance3D(const vector < T1 > & variable1, const vector < T2 > & variable2) {

    double x1 = convert < double > (variable1[0]);
    double y1 = convert < double > (variable1[1]);
    double z1 = convert < double > (variable1[2]);

    double x2 = convert < double > (variable2[0]);
    double y2 = convert < double > (variable2[1]);
    double z2 = convert < double > (variable2[2]);

    return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
}

vector < double > recCircle (const vector < vector < double > > & points, bool qualityCheck = false) {

    double x1 = points[0][0];
    double x2 = points[1][0];
    double x3 = points[2][0];
    double y1 = points[0][1];
    double y2 = points[1][1];
    double y3 = points[2][1];

    double a = x1 - x2;
    double b = y1 - y2;
    double c = x1 - x3;
    double d = y1 - y3;
    double e = ((x1 * x1 - x2 * x2) + (y1 * y1 - y2 * y2)) / 2.0;
    double f = ((x1 * x1 - x3 * x3) + (y1 * y1 - y3 * y3)) / 2.0;

    double det = b * c - a * d;
    if (fabs(det) < 1e-5) {
        if (qualityCheck) return vector < double > {0.0, 0.0, -1.0, 0.0};
        return vector < double > {0.0, 0.0, -1.0};
    }

    double x0 = - (d * e - b * f) / det;
    double y0 = - (a * f - c * e) / det;
    double r0 = distance(vector < double > {x1, y1}, vector < double > {x0, y0});

    if (qualityCheck) {

        double massCenterX = (x1 + x2 + x3) / 3.0;
        double massCenterY = (y1 + y2 + y3) / 3.0;

        double massCenterBias = distance(vector < double > {massCenterX, massCenterY}, vector < double > {x0, y0});
        double quality = massCenterBias / r0;

        return vector < double > {x0, y0, r0, quality};
    }

    return vector < double > {x0, y0, r0};
}

template < typename T1 >
T1 dot(const vector < T1 > & variable1, const vector < T1 > & variable2) {

    auto v1Size = variable1.size();
    if (v1Size != variable2.size()) {
        throw "Error in tooLib::dot()        - - - -         Unequal size of input vectors.\n";
    }

    T1 sum = 0;
    for (int i = 0; i < v1Size; i ++) {
        sum += variable1[i] * variable2[i];
    }

    return sum;
}


// dont use them outside class Array
template < typename T1, typename T2 >
vector < T1 > operator + (const vector < T1 > & variable1, T2 variable2) {

    T1 in = convert < T1 > (variable2);
    auto size = variable1.size();
    vector < T1 > out (size);

    for (int i = 0; i < size; i ++) {
        out[i] = variable1[i] + in;
    }

    return out;
}

template < typename T1, typename T2 >
vector < T1 > operator * (const vector < T1 > & variable1, T2 variable2) {

    T1 in = convert < T1 > (variable2);
    auto size = variable1.size();
    vector < T1 > out (size);

    for (int i = 0; i < size; i ++) {
        out[i] = variable1[i] * in;
    }

    return out;
}

template < typename T1 >
class Array{

private:
    vector < vector < T1 > > arr;
    int row;
    int col;

public:
    Array (int inputRow, int inputCol, T1 val) {

        this -> row = inputRow;
        this -> col = inputCol;
        this -> arr.resize(inputRow);
        for (int iRow = 0; iRow < inputRow; iRow ++) {

            arr[iRow] = vector < T1 > (inputCol, val);
        }
    }

    Array (int inputSize, T1 val) {

        this -> row = inputSize;
        this -> col = inputSize;
        this -> arr.resize(inputSize);
        for (int iRow = 0; iRow < inputSize; iRow ++) {

            arr[iRow] = vector < T1 > (inputSize, 0);
            arr[iRow][iRow] = val;
        }
    }

    Array (const Array < T1 > & inArr) {

        this -> row = inArr.row;
        this -> col = inArr.col;
        this -> arr = inArr.arr;
    }

    vector < int > Size () const { return vector < int > {this -> row, this -> col}; }

    void Show() const {

        for (int iRow = 0; iRow < this -> row; iRow ++) {
            cout << this -> arr[iRow];
        }
    }

    friend ostream & operator << (ostream & os, const Array & inArr) {
        
        inArr.Show();
        return os;
    }

    vector < T1 > & operator [] (unsigned int index) { return arr[index]; }

    const vector < T1 > & operator [] (unsigned int index) const { return arr[index]; }

    vector < T1 > GetRow (unsigned int inputRow) const { return arr[inputRow]; }

    vector < T1 > GetCol (unsigned int inputCol) const {

        vector < T1 > out (this -> row);
        for (int iRow = 0; iRow < this -> row; iRow ++) {
            out[iRow] = arr[iRow][inputCol];
        }

        return out;
    }

    Array < T1 > T () const {

        Array < T1 > out (this -> col, this -> row, 0);

        for (int iRow = 0; iRow < this -> col; iRow ++) {
            out[iRow] = (*this).GetCol(iRow);
        }

        return out;
    }

    void push_back (const vector < T1 > & variable1) {

        this -> row ++;
        this -> arr.push_back(variable1);
    }

    Array < T1 > operator + (const Array < T1 > & variable1) const {


        if (this -> row != variable1.row || this -> col != variable1.col) throw "err\n";
        Array < T1 > out ((*this));

        for (int iRow = 0; iRow < out.row; iRow ++) {
            for (int iCol = 0; iCol < out.col; iCol ++) {
                out[iRow][iCol] += variable1[iRow][iCol];
            }
        }

        return out;
    }

    Array < T1 > operator - (const Array < T1 > & variable1) const { return (*this) + -1.0 * variable1; }

    template < typename T2 >
    Array < T1 > operator + (T2 variable1) const {

        T1 in = convert < T1 > (variable1);
        Array < T1 > out ((*this));

        for (int iRow = 0; iRow < out.row; iRow ++) {
            for (int iCol = 0; iCol < out.col; iCol ++) {
                out[iRow][iCol] += in;
            }
        }

        return out;
    }

    template < typename T2 >
    friend Array < T1 > operator + (T2 variable1, const Array & inArr) { return inArr + variable1; }

    template < typename T2 >
    Array < T1 > operator - (T2 variable1) const { return (*this) + -1 * variable1; }

    template < typename T2 >
    friend Array < T1 > operator - (T2 variable1, const Array & inArr) { return inArr - variable1; }

    template < typename T2 >
    Array < T1 > operator * (T2 variable1) const {

        T1 in = convert < T1 > (variable1);
        Array < T1 > out ((*this));

        for (int iRow = 0; iRow < out.row; iRow ++) {
            for (int iCol = 0; iCol < out.col; iCol ++) {
                out[iRow][iCol] *= in;
            }
        }

        return out;
    }

    template < typename T2 >
    friend Array < T1 > operator * (T2 variable1, const Array & inArr) { return inArr * variable1; }

    friend Array < T1 > Prod (const Array & inArr1, const Array & inArr2) {

        int row1 = inArr1.row;
        int col1 = inArr1.col;
        int row2 = inArr2.row;
        int col2 = inArr2.col;

        if (col1 != row2) return Array < T1 > (0, 0, 0);

        Array < T1 > out (row1, col2, 0);
        for (int iRow = 0; iRow < row1; iRow ++) {
            for (int iCol = 0; iCol < col2; iCol ++) {

                out[iRow][iCol] = dot(inArr1.GetRow(iRow), inArr2.GetCol(iCol));
            }
        }

        return out;
    }
};




namespace RichConst {
    const int RICnrot = 30001; // Rot. matrix no.

    const double RIClgthk_top = 0.02; // LG    top gap
    const double RIClgthk_bot = 0.07; // LG bottom gap


    const double RICotherthk = 0.08; // PMT window thickness
    const double RICcatolength = 1.81; // cathode window length
    const double RICcatogap = 0.03; // Gap btwn PMT pixels
    const double RICeleclength = 2.75; // electronics length below PMT
    const double RICshiheight = 6.5; // This for the new LG with 3.4
    const double RICpmtshield = 0.1; // PMT shield thickness
    const double RICepsln = 0.002; // Epsilon
    const double RICpmtsupportheight = 11.1 - 2.0; // support structure height
    const double PMT_electronics = 3.0; // PMT side size
    const double cato_inner_pixel = 0.42; // Innerr pixel side size in the photocathode
    const double cathode_length = RICcatolength + RICcatogap;
    const double eff_rad_clarity=0.0055;          // clarity used in charge recosntruction
    const int      RICmaxentries = 44;


    const double aeThk                 = 0.1;             // Radiator support thickness (gap btwn radiator tiles).

    const double aglHeight            = 2.5;            // Radiator agl thickness
    const double aglIndex            = 1.0529;        // Index of Agl.
    const double aglSide            = 11.5;            // AGL Radiator tile side length

    const double elecLength         = 2.75;         // Electronics length below PMT.
    const double foilHeight            = 0.1;            // Foil thickness.
    
    const double lgHeight             = 3.0;            // Light guide height without the fixing foil.
    const double lgTopLength        = 3.4;            // Side length of light guide top (Called lg_length in the standalone version)
    const double lgBottomLength     = 1.77;            // Side length on the bottom
    const double lgMirGap             = 0.5;            // Gap length between Mirror and Light Guide.

    const double mirHeight            = 46.32;        // Mirror height.
    const double mirThk             = 0.2;            // Mirror thickness.

    const double NaFHeight            = 0.5;            // NaF radiator thickness
    const double NaFSide            = 8.52;            // NAF Radiator tile side length
    const double NaFIndex             = 1.33;            // Index of NaF.

    const double pmtFoilThk         = 0.1;            // Thickness of the foil over the Light Guide.
    const double pmtLength            = 2.0;            // Phototube length including PMT window.
    const double pmtSupportThk         = 0.6;            // PMT support structure thickness.
    
    const double radHeight            = 3;            // Radiator support structure height.
    const double radMirGap             = 0.1;            // Gap length between Radiator and Mirror.
    const double radPos             = -71.87;        // RICH Radiator position with respect to AMS02.
    const double radRadius            = 60.0;            // RICH Radiator Radius.

    double NaFTransmissionHeight () {
        return 0.5 * NaFHeight + radMirGap + mirHeight + lgMirGap;
    }

    double aglTransmissionHeight () {
        return 0.5 * aglHeight + radMirGap + mirHeight + lgMirGap;
    }

    double pmtTotalHeight() {
        return pmtLength + elecLength + lgHeight + pmtFoilThk;
    }

    double pmtPos() {
        return radHeight + foilHeight + radMirGap + mirHeight +
            lgMirGap + pmtTotalHeight() / 2.0;
    }

}




#endif /* tooLib_h */
