#ifndef FVEC_INCLUDED
#define FVEC_INCLUDED
#include <stdlib.h>

template<typename T, size_t Nm> struct fvec{
	T arr [Nm];
	size_t Ncurr;
	
	fvec():
		Ncurr(0){
	}
	
	void clear(){
		Ncurr = 0;
	}
	
	void push_back(const T& NewEntry){
		arr[Ncurr] = NewEntry;
		Ncurr++;
	}
	
	bool empty() const{
		return (Ncurr == 0);
	}
	
	size_t size() const{
		return Ncurr;
	}
	
	const T& operator[](int ID) const{
		return arr[ID];
	}
	
	T& operator[](int ID){
		return arr[ID];
	}
	
	void erase(int ID){
		Ncurr--;
		arr[ID] = arr[Ncurr];
	}
};

#endif /* FVEC_INCLUDED */
