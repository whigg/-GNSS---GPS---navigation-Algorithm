


#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include<boost/range/numeric.hpp>
#include <iostream>
#include "Header.h"
#include <string>
#include <vector>
#include <algorithm>
#include <memory>
using namespace std;


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