

#define _USE_MATH_DEFINES
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <iostream>
#include <cmath>
#include <string>
#include <vector>

using namespace std;

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

ostream& operator<<(ostream& os, const SatCoord& Scoord) {
	os << "X: " << Scoord.X << " Y: " << Scoord.Y << " Z: " << Scoord.Z;
	return os;
}

struct GpsEphData {
	int    SatNum   = 0;
	int	   TimeOfWeek		= 0;
	int    flag     = 0;
	int    IODC		= 0;
	long   TimeOfColock		= 0;
	int    URA		= 0;
	int    healthS	= 0;
	int	   WeekNum	= 0;
	double TimeGroupDelay		= 0;
	double Af2		= 0;
	double Af1		= 0;
	double Af0		= 0;
	long   TimeOfEphemerids		= 0;
	int	   IODE		= 0; // The Issue of Data, Ephemeris номер набора параметров эфемерид
	double rootBigA	= 0;
	double orbitEccentricity		= 0;
	double M0		= 0;
	double Omega0   = 0;
	double inc0      = 0;
	double argPer	= 0;
	double Deln		= 0;
	double Omegadot = 0;
	double Incdot	= 0;
	double Crc		= 0;
	double Crs		= 0;
	double Cuc		= 0;
	double Cus		= 0;
	double Cic		= 0;
	double Cis		= 0;
	//double cs		= 0;

	GpsEphData() {};
	GpsEphData(const int& _SatNum   ,
			const double& _TimeOfWeek,		
			const double& _flag     ,		
			const double& _IODC		,	
			const double& _TimeOfColock,		
			const double& _URA		,		
			const double& _healthS	,		
			const double& _WeekNum	,
			const double& _TimeGroupDelay,		
			const double& _Af2		,		
			const double& _Af1		,	
			const double& _Af0		,		
			const double& _TimeOfEphemerids,		
			const double& _IODE		,		
			const double& _rootBigA	,		
			const double& _orbitEccentricity,		
			const double& _M0		,		
			const double& _Omega0   ,		
			const double& _inc0     ,	
			const double& _argPer	,	
			const double& _Deln		,		
			const double& _Omegadot , 
			const double& _Incdot	,
			const double& _Crc		,
			const double& _Crs		,
			const double& _Cuc		,
			const double& _Cus		,
			const double& _Cic		,
			const double& _Cis		):      SatNum		(_SatNum	),
											TimeOfWeek			(_TimeOfWeek		),
											flag     	(_flag		),
											IODC		(_IODC		),
											TimeOfColock			(_TimeOfColock		),
											URA			(_URA		),
											healthS		(_healthS	),
											WeekNum		(_WeekNum	),
											TimeGroupDelay			(_TimeGroupDelay		),
											Af2			(_Af2		),
											Af1			(_Af1		),
											Af0			(_Af0		),
											TimeOfEphemerids			(_TimeOfEphemerids		),
											IODE		(_IODE		),
											rootBigA		(_rootBigA		),
											orbitEccentricity			(_orbitEccentricity		),
											M0			(_M0		),
											Omega0   	(_Omega0	),
											inc0     	(_inc0		),
											argPer		(_argPer	),
											Deln		(_Deln		),
											Omegadot 	(_Omegadot	),
											Incdot		(_Incdot	),
											Crc			(_Crc		),
											Crs			(_Crs		),
											Cuc			(_Cuc		),
											Cus			(_Cus		),
											Cic			(_Cic		),
											Cis			(_Cis		)
	
	{};
	void SemiCycle2Radians() {
		M0			*= M_PI;
		Omega0		*= M_PI;
		inc0		*= M_PI;
		argPer		*= M_PI;
		Deln		*= M_PI;
		Omegadot	*= M_PI;
		Incdot		*= M_PI;
	};
};


double getPdelayFromCycle(const double& pdelCycle) {
	return pdelCycle / GPS_L1;
}

SatCoord getSatCoord(const int&	receiverTime, 
					 const double&	pseudoDel, 
					 const GpsEphData& eph,
					 vector<int>& VecDeltTimeGps) {

	// переводим секунды с началад суток в секунды с начала недели
	int weekReceiverTime = receiverTime + 86400 * round((double)(eph.TimeOfWeek - receiverTime) / 86400);

	// вычисляем показания часов спутника на момент предшествия
	int TimeSatPreced = (weekReceiverTime - (int)pseudoDel) % 604800; 
	
	/// Первое приближение, необходимо для вычисления поправки к часам спутника

	// вычисление приращения показаний часов системы GPS на интервале времени
	// между узловым моментом TimeOfEphemerids и моментоп предшествия TimeSatPreced
	int T_nodeIncr = TimeSatPreced - eph.TimeOfEphemerids;

	// вычисление скорректированного среднего движения на момент предшествия TimeSatPreced
	double n = ROOTM / pow(eph.rootBigA, 3) + eph.Deln;

	// Вычисление средней аномалии на момент предшествия
	double MeanAnom = eph.M0 + n*T_nodeIncr;

	/// вычисление поправки к часам спутника
	// вычисление члена коррекции релятивистских эффектов
	double deltTimeRel = CEARTH * eph.orbitEccentricity*eph.rootBigA*sin(MeanAnom);

	// вычисление поправки к часам спутника
	double deltTimeGps = eph.Af0 + eph.Af1*pow((TimeSatPreced - eph.TimeOfColock), 2) + deltTimeRel;

	// сохраняем поправки
	VecDeltTimeGps.push_back(deltTimeGps);
	TimeSatPreced -= deltTimeGps;

	/// Уточнение параметров

	T_nodeIncr = TimeSatPreced - eph.TimeOfEphemerids;
	MeanAnom = eph.M0 + n*T_nodeIncr;

	// Вычисление значения эксцентрической аномалии на момент 
	// предшествия путем решения нелинейного уравнения Кеплера

	double Eccen = MeanAnom;

	double Eccen_next = Eccen - (Eccen - eph.orbitEccentricity*sin(Eccen) - MeanAnom) 
		/ (1 - eph.orbitEccentricity*cos(Eccen));

	while (abs(Eccen_next - Eccen) >= 10e-10) {

		Eccen = Eccen_next;
		Eccen_next = Eccen - (Eccen - eph.orbitEccentricity*sin(Eccen) - MeanAnom) 
			/ (1 - eph.orbitEccentricity*cos(Eccen));
	}

	Eccen = Eccen_next;

	// Вычисление истинной аномалии спутника на момент предшествия
	double TETA = atan2(sqrt(1 - pow(eph.orbitEccentricity, 2))*sin(Eccen), (cos(Eccen) - eph.orbitEccentricity));

	// Вычисление аргумента и поправленного аргумента широты спутника на момент предшествия
	double FI = TETA + eph.argPer;
	double deltaU = eph.Cuc*cos(2 * FI) + eph.Cus*sin(2 * FI);
	double U = FI + deltaU;

	/// Вычисление поправленного радиус вектора спутника на момент предшествия

	// Поправка к радиус вектору
	double deltaRadVec = eph.Crc*cos(2 * FI) + eph.Crs*sin(2 * FI);

	// рвдиус вектор
	double radVec = pow(eph.rootBigA, 2)*(1 - eph.orbitEccentricity*cos(Eccen)) + deltaRadVec;

	/// Вычисление поправленного угла наклонения орбиты на момент предшествия

	// Поправка угла наклонения
	double delta_I_orb = eph.Cic*cos(2 * FI) + eph.Cis*sin(2 * FI);

	// Поправленный угол наклонения
	double I_orbNut = eph.inc0 + delta_I_orb + eph.Incdot*T_nodeIncr;

	// Вычисление поправленной долготы восходящего узла орбиты на момент предшествия
	double Omega = eph.Omega0 + (eph.Omegadot - OMEGA_EARTH)*T_nodeIncr - OMEGA_EARTH*eph.TimeOfEphemerids;

	// Вычисление координат спутника в гринвичской системе координат на момент предшествия
	SatCoord scoord;
	scoord.X = radVec*(cos(U)*cos(Omega) - sin(U)*sin(Omega)*cos(I_orbNut));
	scoord.Y = radVec*(cos(U)*sin(Omega) + sin(U)*cos(Omega)*cos(I_orbNut));
	scoord.Z = radVec*sin(U)*sin(I_orbNut);
	return scoord;
};

pair <vector<SatCoord>, double> SolveLinearEq(const vector<SatCoord>& SatCoords, const vector<SatCoord>& satCoords,
	const double& delay,const vector<int>&  VecDeltTimeGps) {
	boost::numeric::ublas::zero_vector<double> R (SatCoords.size());
	boost::numeric::ublas::zero_vector<double> ResidualVec (SatCoords.size());
	pair <vector<SatCoord>, double> t;
	return t;
};


int main() {
	double receiverTime = 8290;
	auto delay = 130841340.396; // satNum=4
	GpsEphData testGPSEphData{
		 4					 ,
		 180870				 ,
		 20					 ,
		 46					 ,
		 187200				 ,
		 0					 ,
		 0					 ,
		 836				 ,
		 -6.519258e-009		 ,
		 0					 ,
		 -3.865352e-012		 ,
		 -1.586135e-005		 ,
		 187200				 ,
		 46					 ,
		 5153.65454673767	 ,
		 0.0117950872518122	 ,
		 0.164355092216283	 ,
		 0.820717981550843	 ,
		 0.299435918219388	 ,
		 0.347867229022086	 ,
		 1.679382e-009		 ,
		 -2.796583e-009		 ,
		 -1.364242e-010		 ,
		 260.782			 ,
		 -79.03125			 ,
		 -4.248694e-006		 ,
		 5.735084e-006		 ,
		 -1.639128e-007		 ,
		 3.352761e-008
	};
	testGPSEphData.SemiCycle2Radians();
	vector<int> VecDeltTimeGps;
	SatCoord testCoord = getSatCoord(receiverTime,getPdelayFromCycle(delay),testGPSEphData, VecDeltTimeGps);
	cout << testCoord << endl;
	system("pause");
}