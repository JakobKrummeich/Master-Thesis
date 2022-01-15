#ifndef PARTICLE_TYPE_INCLUDED
#define PARTICLE_TYPE_INCLUDED

#include <iostream>
#include <cstdint>

using namespace std;

enum class ParticleType : uint8_t {
	A = 0,
	B = 1
};

ostream& operator<<(ostream& OStream, ParticleType T){
	if (T == ParticleType::A){
		OStream << "A";
	}
	else {
		OStream << "B";
	}
	return OStream;
}

#endif //PARTICLE_TYPE_INCLUDED
