#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#define OMEGA_EARTH 7.2921151467e-5	// angular velocity of earth rotation
#define C			299792458		// speed of light
#define M			3.986005e14		// geocentric gravity const
#define ROOTM		1.996498184322e+7	//
#define	CEARTH		-4.442807633e-10
#define	GPS_L1		1575.42e6		// carrier frequency

struct SatCoord
{
	double X = 0;
	double Y = 0;
	double Z = 0;
};

std::ostream& operator<<(std::ostream& os, const SatCoord& Scoord) {
	os << "X: " << Scoord.X << " Y: " << Scoord.Y << " Z: " << Scoord.Z;
	return os;
}

struct GpsEphData {
	int    SatNum = 0;
	int	   TimeOfWeek = 0;
	int    flag = 0;
	int    IODC = 0;
	long   TimeOfColock = 0;
	int    URA = 0;
	int    healthS = 0;
	int	   WeekNum = 0;
	double TimeGroupDelay = 0;
	double Af2 = 0;
	double Af1 = 0;
	double Af0 = 0;
	long   TimeOfEphemerids = 0;
	int	   IODE = 0; // The Issue of Data, Ephemeris номер набора параметров эфемерид
	double rootBigA = 0;
	double orbitEccentricity = 0;
	double M0 = 0;
	double Omega0 = 0;
	double inc0 = 0;
	double argPer = 0;
	double Deln = 0;
	double Omegadot = 0;
	double Incdot = 0;
	double Crc = 0;
	double Crs = 0;
	double Cuc = 0;
	double Cus = 0;
	double Cic = 0;
	double Cis = 0;
	//double cs		= 0;

	GpsEphData() {};
	GpsEphData(const int& _SatNum,
		const double& _TimeOfWeek,
		const double& _flag,
		const double& _IODC,
		const double& _TimeOfColock,
		const double& _URA,
		const double& _healthS,
		const double& _WeekNum,
		const double& _TimeGroupDelay,
		const double& _Af2,
		const double& _Af1,
		const double& _Af0,
		const double& _TimeOfEphemerids,
		const double& _IODE,
		const double& _rootBigA,
		const double& _orbitEccentricity,
		const double& _M0,
		const double& _Omega0,
		const double& _inc0,
		const double& _argPer,
		const double& _Deln,
		const double& _Omegadot,
		const double& _Incdot,
		const double& _Crc,
		const double& _Crs,
		const double& _Cuc,
		const double& _Cus,
		const double& _Cic,
		const double& _Cis):      SatNum(_SatNum),
		TimeOfWeek(_TimeOfWeek),
		flag(_flag),
		IODC(_IODC),
		TimeOfColock(_TimeOfColock),
		URA(_URA),
		healthS(_healthS),
		WeekNum(_WeekNum),
		TimeGroupDelay(_TimeGroupDelay),
		Af2(_Af2),
		Af1(_Af1),
		Af0(_Af0),
		TimeOfEphemerids(_TimeOfEphemerids),
		IODE(_IODE),
		rootBigA(_rootBigA),
		orbitEccentricity(_orbitEccentricity),
		M0(_M0),
		Omega0(_Omega0),
		inc0(_inc0),
		argPer(_argPer),
		Deln(_Deln),
		Omegadot(_Omegadot),
		Incdot(_Incdot),
		Crc(_Crc),
		Crs(_Crs),
		Cuc(_Cuc),
		Cus(_Cus),
		Cic(_Cic),
		Cis(_Cis)

	{};
	void SemiCycle2Radians() {
		M0 *= M_PI;
		Omega0 *= M_PI;
		inc0 *= M_PI;
		argPer *= M_PI;
		Deln *= M_PI;
		Omegadot *= M_PI;
		Incdot *= M_PI;
	};
};


double getPdelayFromCycle(const double& pdelCycle) {
	return pdelCycle / GPS_L1;
}
