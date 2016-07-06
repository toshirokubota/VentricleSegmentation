#if !defined _myImage_h_
#define _myImage_h_
//#include "stdafx.h"

#include <iostream>
#include <fstream>
#include <cassert>
using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>

#include "mytype.h"

enum pnm_type {IO_PBMA, IO_PGMA, IO_PPMA, IO_PBM, IO_PGM, IO_PPM};

template <class Item>
class myImage {
 public:

  myImage();
  myImage(int bands, int rows, int cols);
  myImage(int bands, int rows, int cols, const Item& init_value);

  ~myImage();

  myImage(const myImage<Item>& image);

  const myImage<Item>& operator =(const myImage<Item>& image);

  // Get the specified pixel value.
  // Each pixel can be indexed by (column, row) or column*row.
  inline Item GetPixel(int band_index, int row_index, int col_index) const {
    assert(Bands>band_index && Rows>row_index && Cols>col_index &&
           band_index>=0 && row_index>=0 && col_index>=0);
    return(pppItem[band_index][row_index][col_index]);
  }
  
  inline Item GetPixel(int pixel_index) const {
    assert(pixel_index>=0 && pixel_index<Bands*Rows*Cols);
    return(pItem[pixel_index]);
  }

  inline Item GetPixelFold(int band_index, int row_index, int col_index) const;
  inline Item GetPixelWrap(int band_index, int row_index, int col_index) const;
  inline Item GetPixelRepeat(int band_index,int row_index,int col_index) const;
  inline Item GetPixelExtrap(int band_index, int row_index, int col_index) const;
  inline Item GetPixelZero(int band_index, int row_index, int col_index) const;
  inline Item GetPixelDefault(int band_index, int row_index, 
                              int col_index, const Item& item) const;

  // Set the specified pixel to the given value.
  inline void SetPixel(int band_index, int row_index, int col_index,\
		const Item& new_pixel) {
    assert(Bands>band_index && Rows>row_index && Cols>col_index &&
           band_index>=0 && row_index>=0 && col_index>=0);
    pppItem[band_index][row_index][col_index] = new_pixel;
  }

  inline void SetPixel(int pixel_index, const Item& new_pixel) {
    assert(pixel_index>=0 && pixel_index<Bands*Rows*Cols);
    pItem[pixel_index] = new_pixel;
  }
  inline void AddPixel(int band_index, int row_index, int col_index,\
		const Item& adder) {
    assert(Bands>band_index && Rows>row_index && Cols>col_index &&
           band_index>=0 && row_index>=0 && col_index>=0);
    pppItem[band_index][row_index][col_index]=
      pppItem[band_index][row_index][col_index]+adder;
  }

  inline void AddPixel(int pixel_index, const Item& adder) {
    assert(pixel_index>=0 && pixel_index<Bands*Rows*Cols);
    pItem[pixel_index] = pItem[pixel_index]+adder;
  }

  // get the protected member data.
  inline int NumPixels() const {return Bands*Rows*Cols;}
  inline int NumBands() const {return Bands;}
  inline int NumRows() const {return Rows;}
  inline int NumCols() const {return Cols;}

  // Artithmetic overloaded operators
  myImage<Item> operator +(const myImage<Item>& image) const;
  void operator +=(const myImage<Item>& image);
  myImage<Item> operator -(const myImage<Item>& image) const;
  void operator -=(const myImage<Item>& image);
  myImage<Item> operator *(const myImage<Item>& image) const;
  void operator *=(const myImage<Item>& image);
  myImage<Item> operator /(const myImage<Item>& image) const;
  void operator /=(const myImage<Item>& image);

  // Logic overloaded operators
  myImage<Item> operator |(const myImage<Item>& image) const;
  void operator |=(const myImage<Item>& image);
  myImage<Item> operator &(const myImage<Item>& image) const;
  void operator &=(const myImage<Item>& image);

  // Scalar arithmetic operators
  myImage<Item> operator +(const Item& value) const;
  void operator +=(const Item& value);
  myImage<Item> operator -(const Item& value) const;
  void operator -=(const Item& value);
  myImage<Item> operator *(const Item& value) const;
  void operator *=(const Item& value);
  myImage<Item> operator ^(const Item& value) const;
  void operator ^=(const Item& value);

  // myImage manipulation operations
  void extractROI(int low_band, int up_band,
                  int low_row, int up_row,
                  int low_col, int up_col);
  myImage<Item> ExtractROI(int low_band, int up_band,
                  int low_row, int up_row,
                  int low_col, int up_col) const;
  void insert(const myImage<Item>& insert,
              int band_offset, int row_offset, int col_offset);
  myImage<Item> Insert(const myImage<Item>& insert,
                     int band_offset, int row_offset, int col_offset);
  void expand(int new_bands, int new_rows, int new_cols,
              int band_offset, int row_offset, int col_offset,
              Item pad_value);
  myImage<Item> Expand(int new_bands, int new_rows, int new_cols,
              int band_offset, int row_offset, int col_offset,
              Item pad_value) const;
  void scaleup(int scale_bands, int scale_rows, int scale_cols);
  myImage<Item> ScaleUp(int scale_bands, int scale_rows, int scale_cols) const;
  void circularShift(int band_shift, int row_shift, int col_shift);
  myImage<Item> CircularShift(int band_shift, int row_shift, int col_shift);

  myImage<unsigned char> Threshold(Item low, Item high, bool inside=true, unsigned char val=1);
  vector<int> Histogram(int numbins, Item low, Item high);
  myImage<unsigned char> Normalize();

	void Clear() {memset((char *)pItem,(char)0,Bands*Rows*Cols*sizeof(Item));}
	void Const(char c) {memset((char *)pItem,(char)c, Bands*Rows*Cols*sizeof(Item));}

	//new version of insert using memcpy
  void insert(myImage<Item>& insert, int offset) {
		assert(Bands*Rows*Cols-offset>=insert.NumPixels());
		memcpy(pItem+offset,insert.GetDataPointer(),insert.NumPixels()*sizeof(Item));
	}
	//dangerous but allows memcpy
	Item* GetDataPointer() const {return pItem;}

  // File IOs
  bool IsPnmFile(const char* filename);
  int ReadPnmFile(const char* filename);
  int WritePnmFile(const char* filename, pnm_type type, int norm) const;
  int ReadRawFile(const char* filename);
  int WriteRawFile(const char* filename) const;

  /*
    Numerical Recipe Interface
    This annoying conversion is the result of Numerical Recepe
    C routines being written originally in Fortran.
    Thus, the array indexing has to be handled with care. I hate it.
  */
  Item** myImage2NRarray(int band) const;
    
  Item* myImage2NRvector(int band, int row) const {
    return pItem+band*Rows*Cols+row*Cols-1;
  }

 protected:
  typedef Item**  ptr2D_Data;
  typedef Item*   ptr1D_Data;

  void ResizeImage(int bands, int rows, int cols);
 private:

  int Bands, Rows, Cols;

  ptr2D_Data* pppItem;
  Item* pItem;

  void _setItemPointers(int bands, int rows, int cols);
  void _resetItemPointers();
};

// Overloading iostream operators for myImage.
template <class Item>
istream& operator>>(istream& s, myImage<Item>& image);

template <class Item>
ostream& operator<<(ostream& s, const myImage<Item>& image);

// Interactive image write routine.
template<class Item>
void
InteractiveImageWrite(const myImage<Item>& image, 
                      char* filename, int normalize);

template<class Item>
void
WriteMultipleBand2PnmFile(const myImage<Item>& image, 
                          char* head, char* extension, int normalize);

template<class Item>
int
CompareImageSize(const myImage<Item>& image1, const myImage<Item>& image2);

template<class Item>
int
ComparePlaneSize(const myImage<Item>& image1, const myImage<Item>& image2);

//typedef myImage<real> RealImage;
//typedef myImage<uchar> ByteImage;
//typedef myImage<int> IntImage;

#include "myImage.cpp"
#include "myImageFileIO.cpp"

#endif /* myImage_h */
